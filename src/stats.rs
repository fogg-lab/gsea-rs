use crate::algorithm::EnrichmentScore;
use crate::utils::Statistic;
use itertools::{izip, Itertools};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GSEASummary {
    pub term: String,
    pub es: f64,
    pub nes: f64,
    pub pval: f64,
    pub fwerp: f64,
    pub fdr: f64,
    pub run_es: Vec<f64>,
    pub hits: Vec<usize>,
    pub esnull: Vec<f64>,
    pub index: Option<usize>,
}

impl GSEASummary {
    fn normalize(&mut self) -> Vec<f64> {
        let e: f64 = self.es;
        let pos_phi: Vec<f64> = self
            .esnull
            .iter()
            .filter_map(|&x| if x >= 0.0 { Some(x) } else { None })
            .collect();
        let neg_phi: Vec<f64> = self
            .esnull
            .iter()
            .filter_map(|&x| if x < 0.0 { Some(x) } else { None })
            .collect();
        let pos_mean = if pos_phi.len() > 0 {
            pos_phi.as_slice().mean()
        } else {
            e
        };
        let neg_mean = if neg_phi.len() > 0 {
            neg_phi.as_slice().mean()
        } else {
            e
        };
        self.nes = if e >= 0.0 {
            e / pos_mean
        } else {
            e / neg_mean.abs()
        };
        let nesnull: Vec<f64> = self
            .esnull
            .iter()
            .map(|&e| {
                if e >= 0.0 {
                    e / pos_mean
                } else {
                    e / neg_mean.abs()
                }
            })
            .collect();
        nesnull
    }

    fn pval(&mut self) {
        let deno: usize;
        let nomi: usize;
        if self.es < 0.0 {
            deno = self.esnull.iter().filter(|&x| *x < 0.0).count();
            nomi = self.esnull.iter().filter(|&x| x <= &self.es).count();
        } else {
            deno = self.esnull.iter().filter(|&x| *x >= 0.0).count();
            nomi = self.esnull.iter().filter(|&x| x >= &self.es).count();
        }
        if deno == 0 {
            self.pval = 1.0;
            return;
        }
        self.pval = (nomi as f64) / (deno as f64);
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GSEAResult {
    pub summaries: Vec<GSEASummary>,
    pub weight: f64,
    pub min_size: usize,
    pub max_size: usize,
    pub nperm: usize,
    nes_concat: Vec<f64>,
    nesnull_concat: Vec<f64>,
    pub seed: u64,
    pub rankings: Vec<Vec<f64>>,
    pub indices: Vec<Vec<usize>>,
}

impl GSEAResult {
    pub fn new(weight: f64, max_size: usize, min_size: usize, nperm: usize, seed: u64) -> Self {
        GSEAResult {
            summaries: Vec::<GSEASummary>::new(),
            weight,
            max_size,
            min_size,
            nperm,
            nes_concat: Vec::<f64>::new(),
            nesnull_concat: Vec::<f64>::new(),
            seed,
            rankings: Vec::<Vec<f64>>::new(),
            indices: Vec::<Vec<usize>>::new(),
        }
    }

    pub fn stat(&mut self, summary: &mut [GSEASummary]) {
        self.nes_concat.clear();
        self.nesnull_concat.clear();
        summary.iter_mut().for_each(|g| {
            g.pval();
            let mut nesnull = g.normalize();
            self.nes_concat.push(g.nes);
            self.nesnull_concat.append(&mut nesnull);
        });
        let fwerps: Vec<f64> = self.fwer_pval();
        let fdrs = self.fdr();
        for (p, q, g) in izip!(fwerps, fdrs, summary) {
            g.fdr = q;
            g.fwerp = p;
        }
        self.nes_concat.clear();
        self.nesnull_concat.clear();
    }

    pub fn _fdr(&mut self) -> Vec<f64> {
        let nes_idx = self.nes_concat.iter().filter(|&x| *x < 0.0).count();
        let fdrs: Vec<f64> = self
            .nes_concat
            .iter()
            .enumerate()
            .map(|(i, &e)| {
                let mut phi_norm: f64;
                let mut phi_obs: f64;
                let mut nes_higher: usize;
                let mut all_higher: usize;
                let mut all_pos: usize;
                let mut nes_pos: usize;
                let mut fdrs_all: Vec<f64> = Vec::new();
                for j in i..self.nperm {
                    let indexes = (j..self.nesnull_concat.len())
                        .step_by(self.nperm)
                        .into_iter();
                    let nesnull: Vec<f64> = indexes.map(|m| self.nesnull_concat[m]).collect();
                    if e < 0.0 {
                        nes_higher = self.nes_concat.iter().filter(|&x| *x <= e).count();
                        all_higher = nesnull.iter().filter(|&x| *x <= e).count();
                        all_pos = nesnull.iter().filter(|&x| *x < 0.0).count();
                        nes_pos = nes_idx;
                    } else {
                        nes_higher = self.nes_concat.iter().filter(|&x| *x >= e).count();
                        all_higher = nesnull.iter().filter(|&x| *x >= e).count();
                        all_pos = nesnull.iter().filter(|&x| *x >= 0.0).count();
                        nes_pos = self.nes_concat.len() - nes_idx;
                    }
                    phi_norm = if all_pos > 0 {
                        (all_higher as f64) / (all_pos as f64)
                    } else {
                        0.0
                    };
                    phi_obs = if nes_pos > 0 {
                        (nes_higher as f64) / (nes_pos as f64)
                    } else {
                        0.0
                    };
                    fdrs_all.push((phi_norm / phi_obs).clamp(f64::MIN, 1.0));
                }
                fdrs_all.as_slice().mean()
            })
            .collect();
        return fdrs;
    }

    pub fn fdr(&mut self) -> Vec<f64> {
        self.nesnull_concat
            .sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        let (indices, nes_sorted) = self.nes_concat.as_slice().argsort(true);
        let all_idx = self.nesnull_concat.partition_point(|x| *x < 0.0);
        let nes_idx = nes_sorted.partition_point(|x| *x < 0.0);
        let fdrs: Vec<f64> = nes_sorted
            .iter()
            .map(|&e| {
                let phi_norm: f64;
                let phi_obs: f64;
                let nes_higher: usize;
                let all_higher: usize;
                let all_pos: usize;
                let nes_pos: usize;
                if e < 0.0 {
                    nes_higher = nes_sorted.partition_point(|x| *x <= e);
                    all_higher = self.nesnull_concat.partition_point(|x| *x <= e);
                    all_pos = all_idx;
                    nes_pos = nes_idx;
                } else {
                    nes_higher = nes_sorted.len() - nes_sorted.partition_point(|x| *x < e);
                    all_higher =
                        self.nesnull_concat.len() - self.nesnull_concat.partition_point(|x| *x < e);
                    all_pos = self.nesnull_concat.len() - all_idx;
                    nes_pos = nes_sorted.len() - nes_idx;
                }
                phi_norm = if all_pos > 0 {
                    (all_higher as f64) / (all_pos as f64)
                } else {
                    0.0
                };
                phi_obs = if nes_pos > 0 {
                    (nes_higher as f64) / (nes_pos as f64)
                } else {
                    0.0
                };
                (phi_norm / phi_obs).clamp(f64::MIN, 1.0)
            })
            .collect();
        let mut fdr_orig_order: Vec<f64> = vec![0.0; fdrs.len()];
        indices.iter().zip(fdrs.iter()).for_each(|(&i, &v)| {
            fdr_orig_order[i] = v;
        });
        return fdr_orig_order;
    }

    fn fwer_pval(&self) -> Vec<f64> {
        let mut max_nes_pos = vec![0.0; self.nperm];
        let mut min_nes_neg = vec![0.0; self.nperm];
        self.nesnull_concat.iter().enumerate().for_each(|(i, &e)| {
            let idx = i % self.nperm;
            if e >= 0.0 {
                max_nes_pos[idx] = e.max(max_nes_pos[idx]);
            } else {
                min_nes_neg[idx] = e.min(min_nes_neg[idx]);
            }
        });
        let fwerp: Vec<f64> = self
            .nes_concat
            .iter()
            .map(|e| {
                if e < &0.0 {
                    (min_nes_neg.iter().filter(|&x| x < e).count() as f64)
                        / (min_nes_neg.iter().filter(|&x| x < &0.0).count() as f64)
                } else {
                    (max_nes_pos.iter().filter(|&x| x >= e).count() as f64)
                        / (max_nes_pos.len() as f64)
                }
            })
            .collect();
        fwerp
    }
}

impl GSEAResult {
    pub fn prerank(&mut self, genes: &[String], metric: &[f64], gmt: &HashMap<&str, &[String]>) {
        let weighted_metric: Vec<f64> = metric.iter().map(|x| x.abs().powf(self.weight)).collect();
        let mut es = EnrichmentScore::new(genes, self.nperm, self.seed);
        let gperm = es.gene_permutation();
        let mut summ = Vec::<GSEASummary>::new();
        for (&term, &gset) in gmt.iter() {
            let gtag = es.gene.isin(gset);
            let gidx = es.hit_index(&gtag);
            if gidx.len() > self.max_size || gidx.len() < self.min_size {
                continue;
            }
            let tag_indicators: Vec<Vec<f64>> = gperm.iter().map(|de| de.isin(&gidx)).collect();
            let (ess, run_es) = es.enrichment_score_gene(&weighted_metric, &tag_indicators);
            let esnull: Vec<f64> = if ess.len() > 1 {
                ess[1..].to_vec()
            } else {
                Vec::new()
            };
            let gss = GSEASummary {
                term: term.to_string(),
                es: ess[0],
                run_es: run_es,
                hits: gidx,
                esnull: esnull,
                ..Default::default()
            };
            summ.push(gss);
        }
        if self.nperm > 0 {
            self.stat(&mut summ);
        }
        self.summaries = summ;
        self.indices.push((0..genes.len()).collect_vec());
        self.rankings.push(metric.to_owned());
    }
}
