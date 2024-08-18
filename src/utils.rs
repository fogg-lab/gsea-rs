use itertools::Itertools;
use rand::seq::SliceRandom;
use rand::Rng;
use std::collections::HashMap;
use std::hash::Hash;

pub trait Statistic {
    fn mean(&self) -> f64;
    fn argsort(&self, ascending: bool) -> (Vec<usize>, Vec<f64>);
}

impl Statistic for &[f64] {
    fn mean(&self) -> f64 {
        let sum = self.iter().sum::<f64>();
        let count = self.len() as f64;
        sum / count
    }
    fn argsort(&self, ascending: bool) -> (Vec<usize>, Vec<f64>) {
        let indices: Vec<usize> = (0..self.len()).collect();
        let sorted_col: Vec<(usize, &f64)> = indices
            .into_iter()
            .zip(self.iter())
            .sorted_by(|&a, &b| a.1.partial_cmp(b.1).unwrap())
            .collect();
        let mut sidx: Vec<usize> = Vec::new();
        let mut sval: Vec<f64> = Vec::new();
        sorted_col.iter().for_each(|(i, &v)| {
            sidx.push(*i);
            sval.push(v);
        });
        if !ascending {
            sidx.reverse();
            sval.reverse();
        }
        (sidx, sval)
    }
}

#[derive(Debug, Clone)]
pub struct DynamicEnum<T> {
    _elt_to_idx: HashMap<T, usize>,
    _idx_to_elt: Vec<T>,
    _num_indices: usize,
}

impl<T> DynamicEnum<T>
where
    T: Eq + Hash + Clone,
{
    pub fn from(vec: &[T]) -> Self {
        let v2m: HashMap<T, usize> = vec
            .iter()
            .enumerate()
            .map(|(i, v)| (v.clone(), i))
            .collect();
        DynamicEnum {
            _num_indices: v2m.len(),
            _elt_to_idx: v2m,
            _idx_to_elt: vec.to_vec(),
        }
    }
    pub fn index_of(&self, element: &T) -> Option<&usize> {
        self._elt_to_idx.get(element)
    }
    pub fn isin(&self, elements: &[T]) -> Vec<f64> {
        let mut _tag_indicator: Vec<f64> = vec![0.0; self._idx_to_elt.len()];
        elements.iter().for_each(|e| {
            if let Some(idx) = self.index_of(e) {
                _tag_indicator[*idx] = 1.0;
            }
        });
        return _tag_indicator;
    }
    pub fn size(&self) -> usize {
        return self._num_indices;
    }
    pub fn shuffle<R>(&mut self, rng: &mut R)
    where
        R: Rng + ?Sized,
    {
        self._idx_to_elt.shuffle(rng);
        self._idx_to_elt.iter().enumerate().for_each(|(i, e)| {
            self._elt_to_idx.insert(e.clone(), i);
        });
    }
}
