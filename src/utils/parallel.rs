#[cfg(feature = "parallel")]
use rayon::prelude::*;

pub(crate) trait MaybeParallelIter<'a, T> {
    type Iter;
    fn maybe_par_iter(&'a self) -> Self::Iter;
}

#[cfg(feature = "parallel")]
impl<'a, T: Sync + 'a> MaybeParallelIter<'a, T> for [T] {
    type Iter = rayon::slice::Iter<'a, T>;

    fn maybe_par_iter(&'a self) -> Self::Iter {
        self.par_iter()
    }
}

#[cfg(not(feature = "parallel"))]
impl<'a, T: 'a> MaybeParallelIter<'a, T> for [T] {
    type Iter = std::slice::Iter<'a, T>;

    fn maybe_par_iter(&'a self) -> Self::Iter {
        self.iter()
    }
}
