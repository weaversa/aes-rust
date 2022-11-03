#[cfg(kani)]
mod aes_kani_proofs {
  use crate::*;

  /**
   * Property demonstrating that Sbox and InvSbox are inverses.
   */
   
  #[kani::proof]
  fn sbox_inverts() {
    let i:u8 = kani::any();
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

  // /**
  //  * Property demonstrating that the step functions of MixColumns and
  //  * InvMixColumns are inverses.
  //  */
   
  // #[kani::proof]
  // fn mixcolumns_step_inverts() {
  //   let s : Word = kani::any();
  //   assert_eq!(inv_mixcolumns_step(mixcolumns_step(s)), s);
  // }

  // /**
  //  * Property demonstrating that MixColumns and InvMixColumns are
  //  * inverses.
  //  */

  // #[kani::proof]
  // fn mixcolumns_inverts() {
  //   let xs : Block = kani::any();
  //   assert_eq!(inv_mixcolumns(mixcolumns(xs)), xs);
  // }

  // /**
  //  * Property demonstrating that Round and InvRound are inverses.
  //  */

  // #[kani::proof]
  // fn round_inverts() {
  //   let xs : Block = kani::any();
  //   let ki : Block = kani::any(); 
  //   assert_eq!(inv_round(round(xs, ki), ki), xs);
  // }

  // /**
  //  * Property demonstrating that Cipher and InvCipher are inverses.
  //  */

  // #[kani::proof]
  // fn cipher_inverts() {
  //   let key       : [Word; 4] = kani::any();
  //   let plaintext : Block     = kani::any();
  //   assert_eq!(inv_cipher(key, cipher(key, plaintext)), plaintext);
  // }
  
}