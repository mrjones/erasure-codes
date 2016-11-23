struct CodeSpec {
    data_chunks: i32,
    local_parities: i32,
    global_parities: i32,
}

type LocalParityIndex = i32;
type DataChunkIndex = i32;
type GlobalParityIndex = i32;

#[derive(Clone, Debug)]
enum ChunkPosition{
    Data(DataChunkIndex),
    LocalParity(LocalParityIndex),
    GlobalParity(GlobalParityIndex),
}

#[derive(Debug)]
pub enum CodeError {
    InvalidCode,
}

impl std::fmt::Display for CodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        match *self {
            CodeError::InvalidCode => {
                try!(write!(f, "Invalid code"));
            },
        }

        return Ok(());
    }
}

impl std::error::Error for CodeError {
    fn description(&self) -> &str {
        match *self {
            CodeError::InvalidCode => "InvalidCode",
        }
    }
}

type CodeResult<T> = std::result::Result<T, CodeError>;

fn local_parity_coverage(spec: &CodeSpec) -> CodeResult<std::collections::HashMap<LocalParityIndex, Vec<DataChunkIndex>>> {
    if spec.data_chunks % spec.local_parities != 0 {
        return Err(CodeError::InvalidCode);
    }

    let data_chunks_per_parity = spec.data_chunks / spec.local_parities;

    let mut result = std::collections::HashMap::new();
    for parity_idx in 0..spec.local_parities {
        let mut v = Vec::new();
        let first_data = parity_idx * data_chunks_per_parity;
        for data_offset in 0..data_chunks_per_parity {
            v.push(first_data + data_offset);
        }
        result.insert(parity_idx, v);
    }
    return Ok(result);
}

fn local_parity_for_data_chunk(spec: &CodeSpec, di: DataChunkIndex) -> DataChunkIndex {
    let data_chunks_per_parity = spec.data_chunks / spec.local_parities;
    return di / data_chunks_per_parity;
}

fn is_recoverable(spec: &CodeSpec, erasures: &Vec<ChunkPosition>) -> bool {
    let needed = spec.data_chunks;

    let mut data_missing = vec![];
    let mut global_parity_missing = vec![];
    let mut unusable_local_parities: Vec<LocalParityIndex> = vec![];
    
    for erasure in erasures {
        match erasure {
            &ChunkPosition::Data(i) => {
                data_missing.push(i);
            },
            &ChunkPosition::GlobalParity(i) => {
                global_parity_missing.push(i);
            },
            &ChunkPosition::LocalParity(i) => {
                unusable_local_parities.push(i);
            },
        }
    }

    let mut locally_recoverable = 0;
    for erasure in erasures {
        match erasure {
            &ChunkPosition::Data(i) => {
                let local_parity_index = local_parity_for_data_chunk(spec, i);
                if !unusable_local_parities.contains(&local_parity_index) {
                    unusable_local_parities.push(local_parity_index);
                    locally_recoverable = locally_recoverable + 1;
                }
            },
            _ => (),
        }
    }

    let data_have = spec.data_chunks - (data_missing.len() as i32);
    let global_parity_have = spec.global_parities - (global_parity_missing.len() as i32);
    return data_have + global_parity_have + locally_recoverable >= needed;
}



fn main() {
    let spec_622 = CodeSpec{
        data_chunks: 6,
        local_parities: 2,
        global_parities: 2,
    };    
}

#[cfg(test)]
mod tests {
    extern crate std;

    use super::ChunkPosition;
    use super::CodeSpec;

    #[test]
    fn test_local_parity_coverage() {
        let spec_622 = CodeSpec {
            data_chunks: 6,
            local_parities: 2,
            global_parities: 2,
        };

        let mut expected = std::collections::HashMap::new();
        expected.insert(0, vec![0, 1, 2]);
        expected.insert(1, vec![3, 4, 5]);
        
        assert_eq!(
            expected,
            super::local_parity_coverage(&spec_622).unwrap());
    }

    #[test]
    fn test_is_recoverable() {
        let spec_622 = CodeSpec {
            data_chunks: 6,
            local_parities: 2,
            global_parities: 2,
        };

        assert!(super::is_recoverable(&spec_622, &vec![]));

        let mut possible_erasures = vec![];
         
        // All 1-failures are recoverable
        for d in 0..spec_622.data_chunks {
            possible_erasures.push(ChunkPosition::Data(d));
        }
        for lp in 0..spec_622.local_parities {
            possible_erasures.push(ChunkPosition::LocalParity(lp));
        }
        for gp in 0..spec_622.global_parities {
            possible_erasures.push(ChunkPosition::GlobalParity(gp));
        }
       
        // All 1-failures are recoverable
        for erasure in &possible_erasures {
            assert!(super::is_recoverable(&spec_622, &vec![erasure.clone()]));
        }

        // All 2-failures are recoverable
        for e1 in 0..possible_erasures.len() {
            for e2 in 0..possible_erasures.len() {
                if e1 != e2 {
                    assert!(super::is_recoverable(
                        &spec_622,
                        &vec![possible_erasures[e1].clone(),
                              possible_erasures[e2].clone()]));
                }
            }
        }

        // All 3-failures are recoverable
        for e1 in 0..possible_erasures.len() {
            for e2 in 0..possible_erasures.len() {
                for e3 in 0..possible_erasures.len() {
                    if e1 != e2 && e1 != e3 && e2 != e3 {
                        assert!(super::is_recoverable(
                            &spec_622,
                            &vec![possible_erasures[e1].clone(),
                                  possible_erasures[e2].clone(),
                                  possible_erasures[e3].clone()]));
                    }
                }
            }
        }        

        // Some 4-failures are recoverable
        assert!(super::is_recoverable(
            &spec_622, &vec![ChunkPosition::Data(0),
                             ChunkPosition::Data(3),
                             ChunkPosition::LocalParity(0),
                             ChunkPosition::LocalParity(1)]));

        assert!(super::is_recoverable(
            &spec_622, &vec![ChunkPosition::Data(0),
                             ChunkPosition::Data(3),
                             ChunkPosition::GlobalParity(0),
                             ChunkPosition::GlobalParity(1)]));

        assert!(super::is_recoverable(
            &spec_622, &vec![ChunkPosition::LocalParity(0),
                             ChunkPosition::LocalParity(1),
                             ChunkPosition::GlobalParity(0),
                             ChunkPosition::GlobalParity(1)]));

        assert!(super::is_recoverable(
            &spec_622, &vec![ChunkPosition::Data(0),
                             ChunkPosition::Data(3),
                             ChunkPosition::LocalParity(0),
                             ChunkPosition::GlobalParity(1)]));

        assert!(super::is_recoverable(
            &spec_622, &vec![ChunkPosition::Data(0),
                             ChunkPosition::Data(3),
                             ChunkPosition::LocalParity(0),
                             ChunkPosition::GlobalParity(0)]));

        // Two erasures in first locality
        assert!(!super::is_recoverable(
            &spec_622, &vec![ChunkPosition::Data(0),
                             ChunkPosition::LocalParity(0),
                             ChunkPosition::GlobalParity(0),
                             ChunkPosition::GlobalParity(1)]));

        assert!(!super::is_recoverable(
            &spec_622, &vec![ChunkPosition::Data(0),
                             ChunkPosition::Data(1),
                             ChunkPosition::GlobalParity(0),
                             ChunkPosition::GlobalParity(1)]));

        // Two erasures in second locality
        assert!(!super::is_recoverable(
            &spec_622, &vec![ChunkPosition::Data(3),
                             ChunkPosition::LocalParity(1),
                             ChunkPosition::GlobalParity(0),
                             ChunkPosition::GlobalParity(1)]));

        assert!(!super::is_recoverable(
            &spec_622, &vec![ChunkPosition::Data(3),
                             ChunkPosition::Data(4),
                             ChunkPosition::GlobalParity(0),
                             ChunkPosition::GlobalParity(1)]));
        /*
        for e1 in 0..possible_erasures.len() {
            for e2 in 0..possible_erasures.len() {
                for e3 in 0..possible_erasures.len() {
                    for e4 in 0..possible_erasures.len() {
                        if e1 != e2 && e1 != e3 && e1 != e4 && e2 != e3 && e2 != e4 && e3 != e4 {

                            
                        assert!(super::is_recoverable(
                            &spec_622,
                            &vec![possible_erasures[e1].clone(),
                                  possible_erasures[e2].clone(),
                                  possible_erasures[e3].clone(),
                                  possible_erasures[e4].clone()]),
                                "Couldn't recover {:?} {:?} {:?} {:?}",
                                possible_erasures[e1].clone(),
                                possible_erasures[e2].clone(),
                                possible_erasures[e3].clone(),
                                possible_erasures[e4].clone());
                        }
                    }
                }
            }
        } 
         */       
    }
}
