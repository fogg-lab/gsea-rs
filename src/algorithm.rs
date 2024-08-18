use crate::utils::DynamicEnum;
use rand::rngs::SmallRng;
use rand::SeedableRng;

pub trait EnrichmentScoreTrait {
    fn running_enrichment_score(&self, metric: &[f64], tag_indicator: &[f64]) -> Vec<f64>;
    fn fast_random_walk(&self, metric: &[f64], tag_indicator: &[f64]) -> f64;
}

#[derive(Debug)]
pub struct EnrichmentScore {
    pub gene: DynamicEnum<String>,
    nperm: usize,
    rng: SmallRng,
}

impl EnrichmentScoreTrait for EnrichmentScore {
    fn running_enrichment_score(&self, metric: &[f64], tag_indicator: &[f64]) -> Vec<f64> {
        let n: f64 = tag_indicator.len() as f64;
        let n_hint: f64 = tag_indicator.iter().sum();
        let n_miss: f64 = n - n_hint;
        let norm_notag: f64 = 1.0 / n_miss;
        let no_tag_indicator: Vec<f64> = tag_indicator.iter().map(|&b| 1.0 - b).collect();
        let sum_correl_tag: Vec<f64> = tag_indicator
            .iter()
            .zip(metric.iter())
            .map(|(&b, &v)| b * v)
            .collect();
        let norm_tag: f64 = 1.0 / sum_correl_tag.iter().sum::<f64>();
        let run_es: Vec<f64> = sum_correl_tag
            .iter()
            .zip(no_tag_indicator.iter())
            .map(|(&b, &v)| b * norm_tag - v * norm_notag)
            .scan(0.0, |acc, x| {
                *acc += x;
                Some(*acc)
            })
            .collect();
        return run_es;
    }

    fn fast_random_walk(&self, metric: &[f64], tag_indicator: &[f64]) -> f64 {
        let ns: f64 = tag_indicator
            .iter()
            .zip(metric.iter())
            .map(|(&b, &v)| b * v)
            .sum::<f64>();
        let n: f64 = metric.len() as f64;
        let k: f64 = tag_indicator.iter().sum::<f64>() as f64;
        let mut res: f64 = 0.0;
        let mut cur: f64 = 0.0;
        let q1: f64 = 1.0 / (n - k);
        let q2: f64 = 1.0 / ns;
        let mut last: f64 = -1.0;
        let p: Vec<f64> = tag_indicator
            .iter()
            .enumerate()
            .filter_map(|(i, &t)| if t > 0.0 { Some(i as f64) } else { None })
            .collect();
        for pos in p {
            cur -= q1 * (pos - last - 1.0);
            if cur.abs() > res.abs() {
                res = cur;
            }
            cur += q2 * metric.get(pos as usize).unwrap();
            if cur.abs() > res.abs() {
                res = cur;
            }
            last = pos;
        }
        return res;
    }
}

impl EnrichmentScore {
    pub fn new(gene: &[String], nperm: usize, seed: u64) -> Self {
        let rng = SmallRng::seed_from_u64(seed);
        EnrichmentScore {
            gene: DynamicEnum::from(gene),
            nperm: nperm + 1,
            rng: rng,
        }
    }

    pub fn hit_index(&self, tag_indicator: &[f64]) -> Vec<usize> {
        tag_indicator
            .iter()
            .enumerate()
            .filter_map(|(i, &t)| if t > 0.0 { Some(i) } else { None })
            .collect()
    }

    pub fn gene_permutation(&mut self) -> Vec<DynamicEnum<usize>> {
        let vec: Vec<usize> = (0..self.gene.size()).collect();
        let mut orig: DynamicEnum<usize> = DynamicEnum::from(&vec);
        let mut gperm: Vec<DynamicEnum<usize>> = Vec::new();
        gperm.push(orig.clone());
        for _ in 1..self.nperm {
            orig.shuffle(&mut self.rng);
            gperm.push(orig.clone());
        }
        return gperm;
    }

    pub fn enrichment_score_gene(
        &mut self,
        metric: &[f64],
        tag_indicators: &[Vec<f64>],
    ) -> (Vec<f64>, Vec<f64>) {
        let es: Vec<f64> = tag_indicators
            .iter()
            .map(|tag| self.fast_random_walk(metric, tag))
            .collect();
        let run_es = self.running_enrichment_score(metric, &tag_indicators[0]);
        return (es, run_es);
    }
}