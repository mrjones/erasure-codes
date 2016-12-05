// TODO: http://web.eecs.utk.edu/~plank/plank/papers/FAST-2013-GF.pdf

// Generalized Pyramid Codes:
// http://research.microsoft.com/en-us/um/people/chengh/papers/pyramid-tos13.pdf
// Page 3:15
pub struct GPCSpec {
    rows: i32,
    cols: i32,
    parities_per_row: i32,
    parities_per_col: i32,
}

pub struct CodeSpec {
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
    MalformedData,
}

impl std::fmt::Display for CodeError {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        match *self {
            CodeError::InvalidCode => {
                try!(write!(f, "Invalid code"));
            },
            CodeError::MalformedData => {
                try!(write!(f, "Malformed data"));
            },
        }

        return Ok(());
    }
}

impl std::error::Error for CodeError {
    fn description(&self) -> &str {
        match *self {
            CodeError::InvalidCode => "InvalidCode",
            CodeError::MalformedData => "MalformedData",
        }
    }
}

type CodeResult<T> = std::result::Result<T, CodeError>;

pub struct CodeWords {
    local_parities: Vec<Vec<u8>>,
    global_parities: Vec<Vec<u8>>,
}

pub fn encode(spec: &CodeSpec, data_chunks: &Vec<Vec<u8>>) -> CodeResult<CodeWords> {
    if spec.data_chunks != (data_chunks.len() as i32) {
        return Err(CodeError::MalformedData);
    }

    let mut code_words = CodeWords{
        local_parities: vec![vec![]; spec.local_parities as usize],
        global_parities: vec![vec![]; spec.global_parities as usize],
    };

    for di in 0..data_chunks.len() {
        let lpi = local_parity_for_data_chunk(spec, di as i32);
        let data_chunk = &data_chunks[di];
        let acc: Vec<u8>;
        if code_words.local_parities[lpi as usize].len() > 0 {
            acc = code_words.local_parities[lpi as usize].clone();
        } else {
            acc = vec![0; data_chunk.len()];
        };

        code_words.local_parities[lpi as usize] =
            acc.iter().zip(data_chunk.iter()).map(|(x1, x2)| x1 ^ x2).collect();
    }

    
    
    return Ok(code_words);
}

// Returns true for entries which will always be zero
fn generate_gzero(spec: &GPCSpec) -> Vec<Vec<bool>> {
    let total_data_chunks = spec.rows * spec.cols;
    let total_row_parities = spec.rows * spec.parities_per_row;
    let total_col_parities = spec.rows * spec.parities_per_col;
    let total_all_chunks = total_data_chunks + total_row_parities + total_col_parities;

    let mut gzero = vec![vec![false; total_data_chunks as usize]; total_all_chunks as usize];
    for row in 0..total_all_chunks {
        for col in 0..total_data_chunks {
            if row < total_data_chunks {
                // Data row
                gzero[row as usize][col as usize] = false;
            } else if row < total_data_chunks + total_row_parities {
                // Row parity
                let row_parity_index = row - total_data_chunks;
                let row_index = row_parity_index / spec.parities_per_row;
                gzero[row as usize][col as usize] =
                    (col < (row_index * spec.cols)) ||
                    (col >= ((row_index + 1) * spec.cols));
            } else {
                // Col parity
                let col_parity_index = row - (total_data_chunks + total_row_parities);
                let col_index = col_parity_index / spec.parities_per_col;
                gzero[row as usize][col as usize] =
                    (col % spec.cols) != col_index;
            }
        }
    }

    return gzero;
}

// TODO(mrjones): u8 is good for GF(2^8). What about bigger fields?
fn global_coefficients(spec: &CodeSpec) -> Vec<Vec<u8>> {
    // Actually Generalized Pyramid Codes:
    // http://research.microsoft.com/en-us/um/people/chengh/papers/pyramid-tos13.pdf
    // Page 3:18
    let k = spec.data_chunks;
    let n = spec.data_chunks + spec.local_parities + spec.global_parities;

    panic!();
}

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
    use super::GPCSpec;

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
    }

    #[test]
    fn test_encode() {
        let spec_622 = CodeSpec {
            data_chunks: 6,
            local_parities: 2,
            global_parities: 2,
        };

        assert!(super::is_recoverable(&spec_622, &vec![]));

        let data = vec![
            vec![0b00000001],
            vec![0b00000010],
            vec![0b00000100],
            vec![0b00001000],
            vec![0b00010000],
            vec![0b00100000],
        ];
        
        let code_word = super::encode(&spec_622, &data).expect("encode");

        assert_eq!(vec![0b00000111], code_word.local_parities[0]);
        assert_eq!(vec![0b00111000], code_word.local_parities[1]);
    }

    #[test]
    fn test_gpc_gzero_2x2_1_1() {
        let spec_2x2_1_1 = GPCSpec{
            rows: 2,
            cols: 2,
            parities_per_row: 1,
            parities_per_col: 1,
        };

        assert_eq!(
            vec![
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, true,  true],
                vec![true,  true,  false, false],
                vec![false, true,  false, true],
                vec![true,  false, true,  false],
            ],
            super::generate_gzero(&spec_2x2_1_1));
    }

    fn test_gpc_gzero_2x2_2_1() {
        let spec_2x2_1_1 = GPCSpec{
            rows: 2,
            cols: 2,
            parities_per_row: 2,
            parities_per_col: 1,
        };

        assert_eq!(
            vec![
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, true,  true],
                vec![false, false, true,  true],
                vec![true,  true,  false, false],
                vec![true,  true,  false, false],
                vec![false, true,  false, true],
                vec![true,  false, true,  false],
            ],
            super::generate_gzero(&spec_2x2_1_1));
    }

    fn test_gpc_gzero_2x2_1_2() {
        let spec_2x2_1_1 = GPCSpec{
            rows: 2,
            cols: 2,
            parities_per_row: 1,
            parities_per_col: 2,
        };

        assert_eq!(
            vec![
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, false, false],
                vec![false, false, true,  true],
                vec![true,  true,  false, false],
                vec![false, true,  false, true],
                vec![false, true,  false, true],
                vec![true,  false, true,  false],
                vec![true,  false, true,  false],
            ],
            super::generate_gzero(&spec_2x2_1_1));
    }
}
