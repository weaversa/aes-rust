#[cfg(test)]
mod aes_tests {
  use crate::*;

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
