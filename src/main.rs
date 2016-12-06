// TODO: http://web.eecs.utk.edu/~plank/plank/papers/FAST-2013-GF.pdf

extern crate la;
extern crate rand;

use rand::Rng;

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

// Should dot be within the field?
// TODO(mrjones): templatize for other fields
fn v_dot(a: &Vec<u8>, b: &Vec<u8>) -> u64 {
    // TODO(mrjones): use zip
    assert_eq!(a.len(), b.len());
    let mut dot = 0;
    for i in 0..a.len() {
        dot += (a[i] as u64) * (b[i] as u64);
    }
    return dot;
}

// TODO(mrjones): u8 overflow?
fn v_add(a: &Vec<u8>, b: &Vec<u8>) -> Vec<u8> {
    return a.iter().zip(b.iter()).map(|(x, y)| x + y).collect();
}

fn v_add_c(a: &Vec<u8>, c: u8) -> Vec<u8> {
    return a.iter().map(|x| x + c).collect();
}

fn v_mult(a: &Vec<u8>, b: &Vec<u8>) -> Vec<u8>{
    return a.iter().zip(b.iter()).map(|(x, y)| x * y).collect();
}

fn v_mult_c(a: &Vec<u8>, c: u8) -> Vec<u8>{
    return a.iter().map(|x| x * c).collect();
}

fn to_matrix(a: &Vec<&Vec<u8>>) -> la::Matrix<f32> {
    assert!(a.len() > 0);

    let mut d = vec![0.0; a.len() * a[0].len()];
    let mut i = 0;
    for r in 0..a.len() {
        for c in 0..a[r].len() {
            d[i] = a[r][c] as f32;
            i = i+1;
        }
    }

    return la::Matrix::new(a.len(), a[0].len(), d);
}

fn rank(svd: &la::SVD<f32>) -> usize {
    let mut rank = 0;
    let count = std::cmp::min(svd.get_s().rows(), svd.get_s().cols());
    for i in 0..count {
        if *svd.get_s().get_ref(i, i) > 0.0 {
            rank = rank + 1
        }
    }

    return rank    
}

fn one_null_basis(svd: &la::SVD<f32>) -> Vec<u8> {
    let count = std::cmp::min(svd.get_s().rows(), svd.get_s().cols());
    assert!(0.1 > *svd.get_s().get_ref(count - 1, count - 1));
    assert!(-0.1 < *svd.get_s().get_ref(count - 1, count - 1));

    let mut null_basis = vec![];
    for r in 0..count {
        null_basis.push(*svd.get_u().get_ref(r, count - 1) as u8);
    }

    return null_basis;
}

fn one_null_basis_wrapper(a: &Vec<&Vec<u8>>) -> Vec<u8> {
    if a.len() == 0 {
        return vec![];
    }
    let m = to_matrix(a);
    let svd = la::SVD::new(&m);
    return one_null_basis(&svd);    
}

fn rank_wrapper(a: &Vec<&Vec<u8>>) -> usize {
    if a.len() == 0 {
        return 0;
    }
    let m = to_matrix(a);

    let svd = la::SVD::new(&m);
    return rank(&svd);
}

fn generate_g(spec: &GPCSpec) -> Vec<Vec<u8>> {
    let total_data_chunks = spec.rows * spec.cols;
    let total_row_parities = spec.rows * spec.parities_per_row;
    let total_col_parities = spec.rows * spec.parities_per_col;
    let total_all_chunks = total_data_chunks + total_row_parities + total_col_parities;

    let mut rng = rand::thread_rng();
    
    // http://research.microsoft.com/en-us/um/people/chengh/papers/pyramid-tos13.pdf
    // Page 3:18
    let mut g = vec![vec![0; total_data_chunks as usize]; total_all_chunks as usize];
    let mut u = vec![vec![0; total_data_chunks as usize]; total_all_chunks as usize];
    let gzero = generate_gzero(spec);

    // G and U start as KxK identity matrix
    for row in 0..total_data_chunks {
        g[row as usize][row as usize] = 1;
        u[row as usize][row as usize] = 1;
    }

    for m in total_data_chunks..total_all_chunks {
        // Randomize row, keeping must-be-zero entries as 0
        for t in 0..total_data_chunks {
            if !gzero[m as usize][t as usize] {
                g[m as usize][t as usize] = rng.gen::<u8>();
            }
        }

        let mut ug = vec![];
        for j in 0..u.len() {
            let d = v_dot(&g[m as usize], &u[j as usize]);
            if d != 0 {
                ug.push(d);
            }

            // TODO(mrjones): what does this mean?
            // if v_dot(u[j], g[m]) === 0, then ug[j] = 0

            let mut ebads = std::collections::HashSet::new();
            let mut uu = vec![];
            for i in 0..j-1 {
                let d = v_dot(&u[i as usize], &u[j as usize]);
                uu.push(d);
                if d != 0 {
                    ebads.insert((ug[i as usize] / uu[i as usize]) as u8);
                }
            }

            let mut e: u8 = 0;
            let e_max = 255 as u8; // TODO(mrjones): compute this based on field size
            let mut found = false;
            for e_candidate in 0..e_max {
                if !ebads.contains(&e_candidate) {
                    e = e_candidate;
                    found = true;
                    break;
                }
            }
            assert!(found);
            
            g[m as usize] = v_add(&g[m as usize], &v_mult_c(&u[j as usize], e));

            for i in 0..j {
                // TODO: What field is this over?
                ug[i] = ug[i] + e as u64 * uu[i];
            }

            // TODO: This is unnecessary, right?
            for t in 0..total_data_chunks {
                if gzero[m as usize][t as usize] {
                    g[m as usize][t as usize] = 0;
                }
            }

            ug[j] = v_dot(&u[j], &g[m as usize]);

            /*
            for si in 0..g.len() {
                for sj in si..g.len() {
                    // S' is G withour row si or sj:

                    // TODO(mrjones): don't copy so much!
                    let s_prime = vec![];
                    for i in 0..g.len() {
                        if i != si && i != sj {
                            s_prime.push(&g[i]);
                        }
                    }
                    s_prime.push(&g[m as usize]);

                    if rank(&s_prime) == (total_data_chunks - 1) as usize {
                        u.push(null_space(s_prime));
                    }
                }
            }
             */
        }
    }


    return g;
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

    
    #[test]
    fn test_rank_1() {
        let r1 = vec![1, 0];
        let r2 = vec![0, 1];
        let m = vec![&r1, &r2];

        assert_eq!(2, super::rank_wrapper(&m));
    }

    #[test]
    fn test_rank_2() {
        let r1 = vec![1, 0];
        let r2 = vec![1, 0];
        let m = vec![&r1, &r2];

        assert_eq!(1, super::rank_wrapper(&m));
    }

    #[test]
    fn test_rank_3() {
        let r1 = vec![1, 0, 0];
        let r2 = vec![0, 1, 0];
        let r3 = vec![0, 2, 0];
        let m = vec![&r1, &r2, &r3];

        assert_eq!(2, super::rank_wrapper(&m));
    }

    #[test]
    fn test_one_null_basis() {
        let r1 = vec![1, 0, 0];
        let r2 = vec![0, 1, 0];
        let r3 = vec![0, 0, 0];
        let m = vec![&r1, &r2, &r3];

        assert_eq!(vec![0, 0, 1], super::one_null_basis_wrapper(&m));
    }
}
