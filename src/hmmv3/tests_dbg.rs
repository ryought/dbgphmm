#[cfg(test)]
mod tests {
    use super::*;
    use crate::common::{ei, ni, Reads};
    use crate::dbg::mocks::*;
    use crate::hmmv2::mocks::mock_linear_phmm;
    use crate::hmmv2::params::PHMMParams;
    use crate::io::cytoscape::NodeAttrVec;
    use crate::prob::p;
    use crate::vector::DenseStorage;

    #[test]
    fn hmm_dbg_changing_edge_copy_nums() {
        let mut param = PHMMParams::default();
        let reads = Reads {
            reads: vec![b"AACTAGCTT".to_vec()],
        };

        let mut dbg = mock_intersection();
        let phmm = dbg.to_phmm(PHMMParams::default());
        let nf = phmm.to_node_freqs(&reads);
        println!("{}", nf);
        let e = 0.9;
        assert!(nf[ni(1)] > e);
        assert!(nf[ni(2)] > e);
        assert!(nf[ni(8)] > e);
        assert!(nf[ni(12)] > e);
        assert!(nf[ni(13)] > e);
        assert!(nf[ni(16)] > e);
        assert!(nf[ni(20)] > e);
        assert!(nf[ni(21)] > e);

        // phmm.draw_node_vec(&nf);
        // println!(
        //     "{}",
        //     dbg.to_cytoscape_with_attrs(&[NodeAttrVec::Freq(nf)], &[])
        // );

        // let (ncn, ecn) = dbg.to_copy_nums_of_seq(b"AACTAGCTT").unwrap();
        let (ncn, ecn) = dbg.to_copy_nums_of_seq(b"CCGTAGGGC").unwrap();
        dbg.set_node_copy_nums(&ncn);
        dbg.set_edge_copy_nums(Some(&ecn));
        // println!("{}", dbg);
        let phmm = dbg.to_phmm(PHMMParams::default());
        let nf = phmm.to_node_freqs(&reads);
        println!("{}", nf);
        let e = 0.001;
        assert!(nf[ni(1)] < e);
        assert!(nf[ni(2)] < e);
        assert!(nf[ni(8)] < e);
        assert!(nf[ni(12)] < e);
        assert!(nf[ni(13)] < e);
        assert!(nf[ni(16)] < e);
        assert!(nf[ni(20)] < e);
        assert!(nf[ni(21)] < e);
    }
}
