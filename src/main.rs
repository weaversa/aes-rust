/**
 * This sequence of bytes comprises the input, output, State, and
 * Round Key [FIPS-PUB-197], Section 2.1.
 */

type Block = [u8; 16];

/**
 * This sequence of bytes represents a 32-bit word [FIPS-PUB-197],
 * Section 2.1.
 */

type Word = [u8; 4];

/**
 * The round constant word array [FIPS-PUB-197], Section 2.1 and
 * 5.2. Constants are use here rather than computing the values in
 * place.
 */

static RCON : [Word; 11] =
       [ [0x00, 0x00, 0x00, 0x00]  // undefined
       , [0x01, 0x00, 0x00, 0x00], [0x02, 0x00, 0x00, 0x00]
       , [0x04, 0x00, 0x00, 0x00], [0x08, 0x00, 0x00, 0x00]
       , [0x10, 0x00, 0x00, 0x00], [0x20, 0x00, 0x00, 0x00]
       , [0x40, 0x00, 0x00, 0x00], [0x80, 0x00, 0x00, 0x00]
       , [0x1b, 0x00, 0x00, 0x00], [0x36, 0x00, 0x00, 0x00]
       ];

/**
 * Multiplication by x (i.e., 0b00000010 or 0x02) can be implemented
 * at the byte level as a left shift and a subsequent conditional
 * bitwise XOR with 0x1b [FIPS-PUB-197], Section 4.2.1.
 */

fn xtime (x : u8) -> u8 {
    x << 1 ^ (if x > 127 { 0x1b } else { 0x00 })
}

/**
 * Multiplication of two bytes representing coefficients of
 * polynomials in GF 2^^8 [FIPS-PUB-197], Section 4.3.
 */

fn mul_bytes (x : u8, y : u8) -> u8 {
    let mut a : u8 = 0x00;
    for i in 0..=7 {
      a = xtime(a) ^ (if x>>(7-i) & 0x01 == 0x01 { y } else { 0x00 });
    }
    return a;
}

/**
 * S-box [FIPS-PUB-197], Section 5.1.1, Figure 7.
 */

static SBOX : [u8; 256] =
    [ 0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76
    , 0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0
    , 0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15
    , 0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75
    , 0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84
    , 0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf
    , 0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8
    , 0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2
    , 0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73
    , 0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb
    , 0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79
    , 0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08
    , 0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a
    , 0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e
    , 0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf
    , 0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16 ];

/**
 * Inverse S-box [FIPS-PUB-197], Section 5.3.2, Figure 14.
 */

static INVSBOX : [u8; 256] =
    [ 0x52, 0x09, 0x6a, 0xd5, 0x30, 0x36, 0xa5, 0x38, 0xbf, 0x40, 0xa3, 0x9e, 0x81, 0xf3, 0xd7, 0xfb
    , 0x7c, 0xe3, 0x39, 0x82, 0x9b, 0x2f, 0xff, 0x87, 0x34, 0x8e, 0x43, 0x44, 0xc4, 0xde, 0xe9, 0xcb
    , 0x54, 0x7b, 0x94, 0x32, 0xa6, 0xc2, 0x23, 0x3d, 0xee, 0x4c, 0x95, 0x0b, 0x42, 0xfa, 0xc3, 0x4e
    , 0x08, 0x2e, 0xa1, 0x66, 0x28, 0xd9, 0x24, 0xb2, 0x76, 0x5b, 0xa2, 0x49, 0x6d, 0x8b, 0xd1, 0x25
    , 0x72, 0xf8, 0xf6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xd4, 0xa4, 0x5c, 0xcc, 0x5d, 0x65, 0xb6, 0x92
    , 0x6c, 0x70, 0x48, 0x50, 0xfd, 0xed, 0xb9, 0xda, 0x5e, 0x15, 0x46, 0x57, 0xa7, 0x8d, 0x9d, 0x84
    , 0x90, 0xd8, 0xab, 0x00, 0x8c, 0xbc, 0xd3, 0x0a, 0xf7, 0xe4, 0x58, 0x05, 0xb8, 0xb3, 0x45, 0x06
    , 0xd0, 0x2c, 0x1e, 0x8f, 0xca, 0x3f, 0x0f, 0x02, 0xc1, 0xaf, 0xbd, 0x03, 0x01, 0x13, 0x8a, 0x6b
    , 0x3a, 0x91, 0x11, 0x41, 0x4f, 0x67, 0xdc, 0xea, 0x97, 0xf2, 0xcf, 0xce, 0xf0, 0xb4, 0xe6, 0x73
    , 0x96, 0xac, 0x74, 0x22, 0xe7, 0xad, 0x35, 0x85, 0xe2, 0xf9, 0x37, 0xe8, 0x1c, 0x75, 0xdf, 0x6e
    , 0x47, 0xf1, 0x1a, 0x71, 0x1d, 0x29, 0xc5, 0x89, 0x6f, 0xb7, 0x62, 0x0e, 0xaa, 0x18, 0xbe, 0x1b
    , 0xfc, 0x56, 0x3e, 0x4b, 0xc6, 0xd2, 0x79, 0x20, 0x9a, 0xdb, 0xc0, 0xfe, 0x78, 0xcd, 0x5a, 0xf4
    , 0x1f, 0xdd, 0xa8, 0x33, 0x88, 0x07, 0xc7, 0x31, 0xb1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xec, 0x5f
    , 0x60, 0x51, 0x7f, 0xa9, 0x19, 0xb5, 0x4a, 0x0d, 0x2d, 0xe5, 0x7a, 0x9f, 0x93, 0xc9, 0x9c, 0xef
    , 0xa0, 0xe0, 0x3b, 0x4d, 0xae, 0x2a, 0xf5, 0xb0, 0xc8, 0xeb, 0xbb, 0x3c, 0x83, 0x53, 0x99, 0x61
    , 0x17, 0x2b, 0x04, 0x7e, 0xba, 0x77, 0xd6, 0x26, 0xe1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0c, 0x7d ];

/**
 * SubBytes transformation [FIPS-PUB-197], Section 5.1.1.
 */

fn subbytes (xs : Block) -> Block {
  let mut ret : Block = [0; 16];
  for i in 0..=15 {
    ret[i] = SBOX[xs[i] as usize];
  }
  return ret;
}

/**
 * InvSubBytes transformation [FIPS-PUB-197], Section 5.3.2.
 */

fn inv_subbytes (xs : Block) -> Block {
  let mut ret : Block = [0; 16];
  for i in 0..=15 {
    ret[i] = INVSBOX[xs[i] as usize];
  }
  return ret;
}

/**
 * ShiftRows transformation [FIPS-PUB-197], Section 5.1.2.
 */

fn shiftrows (xs : Block) -> Block {
  static PERM : [usize; 16] = [0, 5, 10, 15, 4, 9, 14, 3, 8, 13, 2, 7, 12, 1, 6, 11];
  let mut ret : Block = [0; 16];
  for i in 0..=15 {
    ret[i] = xs[PERM[i]];
  }
  return ret;
}

/**
 * InvShiftRows transformation [FIPS-PUB-197], Section 5.3.1.
 */

fn inv_shiftrows (xs : Block) -> Block {
  static PERM : [usize; 16] = [0, 13, 10, 7, 4, 1, 14, 11, 8, 5, 2, 15, 12, 9, 6, 3];
  let mut ret : Block = [0; 16];
  for i in 0..=15 {
    ret[i] = xs[PERM[i]];
  }
  return ret;
}

/**
 * One step of MixColumns [FIPS-PUB-197], Section 5.1.3.
 */
 
fn mixcolumns_step ([s0, s1, s2, s3] : Word) -> Word {
  let t0 = mul_bytes(0x02, s0) ^ mul_bytes(0x03, s1) ^ s2                  ^ s3;
  let t1 = s0                  ^ mul_bytes(0x02, s1) ^ mul_bytes(0x03, s2) ^ s3;
  let t2 = s0                  ^ s1                  ^ mul_bytes(0x02, s2) ^ mul_bytes(0x03, s3);
  let t3 = mul_bytes(0x03, s0) ^ s1                  ^ s2                  ^ mul_bytes(0x02, s3);
  return [t0, t1, t2, t3];
}

/**
 * MixColumns transformation [FIPS-PUB-197], Section 5.1.3.
 */

fn mixcolumns ([s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15] : Block) -> Block {
  let [t0 , t1 , t2 , t3 ] = mixcolumns_step([s0 , s1 , s2 , s3 ]);
  let [t4 , t5 , t6 , t7 ] = mixcolumns_step([s4 , s5 , s6 , s7 ]);
  let [t8 , t9 , t10, t11] = mixcolumns_step([s8 , s9 , s10, s11]);
  let [t12, t13, t14, t15] = mixcolumns_step([s12, s13, s14, s15]);
  return [t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15];
}

/**
 * One step of InvMixColumns [FIPS-PUB-197], Section 5.3.3.
 */

fn inv_mixcolumns_step ([s0, s1, s2, s3] : Word) -> Word {
  let t0 = mul_bytes(0x0e, s0) ^ mul_bytes(0x0b, s1) ^ mul_bytes(0x0d, s2) ^ mul_bytes(0x09, s3);
  let t1 = mul_bytes(0x09, s0) ^ mul_bytes(0x0e, s1) ^ mul_bytes(0x0b, s2) ^ mul_bytes(0x0d, s3);
  let t2 = mul_bytes(0x0d, s0) ^ mul_bytes(0x09, s1) ^ mul_bytes(0x0e, s2) ^ mul_bytes(0x0b, s3);
  let t3 = mul_bytes(0x0b, s0) ^ mul_bytes(0x0d, s1) ^ mul_bytes(0x09, s2) ^ mul_bytes(0x0e, s3);
  return [t0, t1, t2, t3];
}

/**
 * InvMixColumns transformation [FIPS-PUB-197], Section 5.3.3.
 */

fn inv_mixcolumns ([s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15] : Block) -> Block {
  let [t0 , t1 , t2 , t3 ] = inv_mixcolumns_step([s0 , s1 , s2 , s3 ]);
  let [t4 , t5 , t6 , t7 ] = inv_mixcolumns_step([s4 , s5 , s6 , s7 ]);
  let [t8 , t9 , t10, t11] = inv_mixcolumns_step([s8 , s9 , s10, s11]);
  let [t12, t13, t14, t15] = inv_mixcolumns_step([s12, s13, s14, s15]);
  return [t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14, t15];
}

/**
 * Performs one round of AES encryption [FIPS-PUB-197], Section 5.1.
 */

fn round (mut xs : Block, ki : Block) -> Block {
  xs = mixcolumns(shiftrows(subbytes(xs)));
  for i in 0..=15 {
    xs[i] = xs[i] ^ ki[i];
  }
  return xs;
}

/**
 * Performs one round of AES decryption [FIPS-PUB-197], Section 5.3.
 */

fn inv_round (mut xs : Block, ki : Block) -> Block {
  for i in 0..=15 {
    xs[i] = xs[i] ^ ki[i];
  }
  return inv_subbytes(inv_shiftrows(inv_mixcolumns(xs)));
}

/**
 * Takes a four-byte input word and applies the S-box to each of the
 * four bytes to produce an output word [FIPS-PUB-197], Section 5.2.
 */

fn subword (w : Word) -> Word {
  let mut ret = [0, 0, 0, 0];
  for i in 0..=3 {
     ret[i] = SBOX[w[i] as usize];
  }
  return ret;
}

/**
 * Takes a word [a0,a1,a2,a3] as input, performs a cyclic permutation,
 * and returns the word [a1,a2,a3,a0] [FIPS-PUB-197], Section 5.2.
 */

fn rotword ([a0, a1, a2, a3] : Word) -> Word {
  [a1, a2, a3, a0]
}

fn xorw (mut a : Word, b : Word) -> Word {
  for i in 0..=3 {
    a[i] ^= b[i]
  }
  return a;
}

/**
 * The main AES key expansion routine [FIPS-PUB-197], Section 5.2.
 */

fn key_expansion (k : [Word; 4]) -> [Block; 11] {
  let nk = 128 / 32;
  let nr = nk + 6;
  let mut w = [[0; 4]; 44];
  for i in 0..=4 * (nr + 1) - 1 {
    if i < nk { w[i] = k[i] }
    else if i % nk == 0 { w[i] = xorw(w[i - nk], xorw(subword(rotword(w[i-1])), RCON[i / nk])) }
    else if (i % nk == 4) & (nk > 6) { w[i] = xorw(w[i-nk], subword(w[i-1])) }
    else { w[i] = xorw(w[i-nk], w[i-1]) };
  }
  let mut bs = [[0u8; 16]; 11];
  for i in 0..=10 {
    for j in 0..=15 {
      bs[i][j] = w[i * 4 + j / 4][j%4];
    }
  }
  return bs;
}

fn xorb (mut a : Block, b : Block) -> Block {
  for i in 0..=15 {
    a[i] ^= b[i]
  }
  return a;
}

/**
 * The Cipher function performs AES encryption [FIPS-PUB-197], Section
 * 5.1.
 */

pub fn cipher (key : [Word; 4], plaintext : Block) -> Block {
  let ks  = key_expansion(key);
  let pre = xorb(plaintext, ks[0]);
  let mut mid = pre;
  for i in 1..=9 {
    mid = round(mid, ks[i]);
  }
  let ciphertext = xorb(shiftrows(subbytes(mid)), ks[10]);
  return ciphertext;
}

/**
 * The InvCipher function performs AES decryption [FIPS-PUB-197],
 * Section 5.1.
 */

pub fn inv_cipher (key : [Word; 4], ciphertext : Block) -> Block {
  let ks  = key_expansion(key);
  let pre = inv_subbytes(inv_shiftrows(xorb(ciphertext, ks[10])));
  let mut mid = pre;
  for i in 1..=9 {
    mid = inv_round(mid, ks[10-i]);
  }
  let plaintext = xorb(mid, ks[0]);
  return plaintext;
}

#[cfg(test)]
mod aes_tests {
  use super::*;

  /**
   * Property demonstrating that Sbox and InvSbox are inverses.
   */
   
  #[test]
  fn sbox_inverts() {
    for i in 0..=255 {
      assert_eq!(INVSBOX[SBOX[i as usize] as usize], i);
    }
  }

  /**
   * Property demonstrating that SubBytes and InvSubBytes are inverses.
   */
   
  #[test]
  fn subbytes_inverts() {
    for _i in 0..=127 {  // Run 128 random tests of 2^^128 possible tests
      let xs : Block = rand::random(); 
      assert_eq!(inv_subbytes(subbytes(xs)), xs);
    }
  }

  /**
   * Property demonstrating that ShiftRows and InvShiftRows are
   * inverses.
   */
   
  #[test]
  fn shiftrows_inverts() {
    for _i in 0..=127 {  // Run 128 random tests of 2^^128 possible tests
      let xs : Block = rand::random(); 
      assert_eq!(inv_shiftrows(shiftrows(xs)), xs);
    }
  }

  /**
   * Property demonstrating that the step functions of MixColumns and
   * InvMixColumns are inverses.
   */
   
  #[test]
  fn mixcolumns_step_inverts() {
    for _i in 0..=127 {  // Run 128 random tests of 2^^32 possible tests
      let s : Word = rand::random(); 
      assert_eq!(inv_mixcolumns_step(mixcolumns_step(s)), s);
    }
  }

  /**
   * Property demonstrating that MixColumns and InvMixColumns are
   * inverses.
   */

  #[test]
  fn mixcolumns_inverts() {
    for _i in 0..=127 {  // Run 128 random tests of 2^^128 possible tests
      let xs : Block = rand::random(); 
      assert_eq!(inv_mixcolumns(mixcolumns(xs)), xs);
    }
  }

  /**
   * Property demonstrating that Round and InvRound are inverses.
   */

  #[test]
  fn round_inverts() {
    for _i in 0..=127 {  // Run 128 random tests of 2^^256 possible tests
      let xs : Block = rand::random();
      let ki : Block = rand::random(); 
      assert_eq!(inv_round(round(xs, ki), ki), xs);
    }
  }

  /**
   * Property demonstrating that Cipher and InvCipher are inverses.
   */

  #[test]
  fn cipher_inverts() {
    for _i in 0..=127 {  // Run 128 random tests of 2^^256 possible tests
      let key   : [Word; 4] = rand::random();
      let plaintext : Block = rand::random(); 
      assert_eq!(inv_cipher(key, cipher(key, plaintext)), plaintext);
    }
  }

  #[test]
  fn test_vectors() {
    let plaintext : Block  = [ 0x00, 0x11, 0x22, 0x33
                             , 0x44, 0x55, 0x66, 0x77
                             , 0x88, 0x99, 0xaa, 0xbb
                             , 0xcc, 0xdd, 0xee, 0xff ];
    let key : [Word; 4]    = [ [0x00, 0x01, 0x02, 0x03]
                             , [0x04, 0x05, 0x06, 0x07]
                             , [0x08, 0x09, 0x0a, 0x0b]
                             , [0x0c, 0x0d, 0x0e, 0x0f] ];
    let ciphertext : Block = [ 0x69, 0xc4, 0xe0, 0xd8
                             , 0x6a, 0x7b, 0x04, 0x30
                             , 0xd8, 0xcd, 0xb7, 0x80
                             , 0x70, 0xb4, 0xc5, 0x5a ];
    assert_eq!(cipher(key, plaintext), ciphertext);
    assert_eq!(inv_cipher(key, ciphertext), plaintext);
  }
}

#[cfg(kani)]
mod aes_kani_proofs {
  use super::*;

  /**
   * Property demonstrating that Sbox and InvSbox are inverses.
   */
   
  #[kani::proof]
  fn sbox_inverts() {
    let i = kani::any();
    
    kani::assume(i <= 255);
    /* The above line gives the warning:
     *   > warning: comparison is useless due to type limits
     * However, leaving it out gives an error:
     *   > error: error: conflicting types for variable
     *     'sbox_inverts::1::var_12'
     */
     
    assert_eq!(INVSBOX[SBOX[i as usize] as usize], i);
  }

  /**
   * Property demonstrating that SubBytes and InvSubBytes are inverses.
   */
   
  #[kani::proof]
  fn subbytes_inverts() {
    let xs : Block = kani::any();
    assert_eq!(inv_subbytes(subbytes(xs)), xs);
  }

  /**
   * Property demonstrating that ShiftRows and InvShiftRows are
   * inverses.
   */
   
  #[kani::proof]
  fn shiftrows_inverts() {
    let xs : Block = kani::any();
    assert_eq!(inv_shiftrows(shiftrows(xs)), xs);
  }

  /**
   * Property demonstrating that the step functions of MixColumns and
   * InvMixColumns are inverses.
   */
   
  #[kani::proof]
  fn mixcolumns_step_inverts() {
    let s : Word = kani::any();
    assert_eq!(inv_mixcolumns_step(mixcolumns_step(s)), s);
  }

  /**
   * Property demonstrating that MixColumns and InvMixColumns are
   * inverses.
   */

  #[kani::proof]
  fn mixcolumns_inverts() {
    let xs : Block = kani::any();
    assert_eq!(inv_mixcolumns(mixcolumns(xs)), xs);
  }

  /**
   * Property demonstrating that Round and InvRound are inverses.
   */

  #[kani::proof]
  fn round_inverts() {
    let xs : Block = kani::any();
    let ki : Block = kani::any(); 
    assert_eq!(inv_round(round(xs, ki), ki), xs);
  }

  /**
   * Property demonstrating that Cipher and InvCipher are inverses.
   */

  #[kani::proof]
  fn cipher_inverts() {
    let key       : [Word; 4] = kani::any();
    let plaintext : Block     = kani::any();
    assert_eq!(inv_cipher(key, cipher(key, plaintext)), plaintext);
  }

}

#[cfg(crux)]
mod aes_crux_proofs {
  extern crate crucible;
  use super::*;
  use crucible::*;

  /**
   * Property demonstrating that Sbox and InvSbox are inverses.
   */
   
  #[crux_test]
  fn sbox_inverts() {
    let i = u8::symbolic("i");
    assert_eq!(INVSBOX[SBOX[i as usize] as usize], i);
  }

  /**
   * Property demonstrating that SubBytes and InvSubBytes are inverses.
   */
   
  #[crux_test]
  fn subbytes_inverts() {
    let mut xs : Block = [0; 16];
    for i in 0..=15 {
      let j = u8::symbolic("xs_i");
      crucible_assume!(INVSBOX[SBOX[j as usize] as usize] == j);
      xs[i] = j;
    }
    assert_eq!(inv_subbytes(subbytes(xs)), xs);
  }

  /**
   * Property demonstrating that ShiftRows and InvShiftRows are
   * inverses.
   */
   
  #[crux_test]
  fn shiftrows_inverts() {
    let mut xs : Block = [0; 16];
    for i in 0..=15 {
      xs[i] = u8::symbolic("xs_i");
    }
    assert_eq!(inv_shiftrows(shiftrows(xs)), xs);
  }

  /**
   * Property demonstrating that the step functions of MixColumns and
   * InvMixColumns are inverses.
   */
   
  #[crux_test]
  fn mixcolumns_step_inverts() {
    let mut s : Word = [0; 4];
    for i in 0..=3 {
      s[i] = u8::symbolic("s_i");
    }
    assert_eq!(inv_mixcolumns_step(mixcolumns_step(s)), s);
  }

  /**
   * Property demonstrating that MixColumns and InvMixColumns are
   * inverses.
   */

  #[crux_test]
  fn mixcolumns_inverts() {
    let mut xs : Block = [0; 16];
    for i in 0..=15 {
      xs[i] = u8::symbolic("xs_i");
    }
    let [s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, s13, s14, s15] : Block = xs;
    crucible_assume!(
      inv_mixcolumns_step(mixcolumns_step([s0 , s1 , s2 , s3 ]))
                                       == [s0 , s1 , s2 , s3 ]);
    crucible_assume!(
      inv_mixcolumns_step(mixcolumns_step([s4 , s5 , s6 , s7 ]))
                                       == [s4 , s5 , s6 , s7 ]);
    crucible_assume!(
      inv_mixcolumns_step(mixcolumns_step([s8 , s9 , s10, s11]))
                                       == [s8 , s9 , s10, s11]);
    crucible_assume!(
      inv_mixcolumns_step(mixcolumns_step([s12, s13, s14, s15]))
                                       == [s12, s13, s14, s15]);
    assert_eq!(inv_mixcolumns(mixcolumns(xs)), xs);
  }

  /**
   * Property demonstrating that Round and InvRound are inverses.
   */

  #[crux_test]
  fn round_inverts() {
    let mut xs : Block = [0; 16];
    let mut ki : Block = [0; 16];
    for i in 0..=15 {
      xs[i] = u8::symbolic("xs_i");
      ki[i] = u8::symbolic("ki_i");
    }

    assert_eq!(inv_round(round(xs, ki), ki), xs);
  }

}

fn main() {
  for i in 0..=0xffffffffu64 {
    let s : Word = [ (i     & 0xff) as u8
                   , (i>> 8 & 0xff) as u8
                   , (i>>16 & 0xff) as u8
                   , (i>>24 & 0xff) as u8];
    assert_eq!(inv_mixcolumns_step(mixcolumns_step(s)), s);
    if i % 1000000 == 0 { println!("i = {}", i) };
  }
}
