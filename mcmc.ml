(* 

   mcmc.ml
   Markov Chain Monte Carlo

   Andrei de A. Formiga, 2013-03-22

*)

type chain = {
    lud: Matrix.t -> float;     (* log unnormalized density: given state, gives density *)
    mutable state: Matrix.t;    (* current state of the chain *)
    f: int -> float * Matrix.t; (* simulation function: given the number of iterations, return accept. rate, result matrix *)
    outf: Matrix.t -> Matrix.t
}


let id x = x

(** 
    Metropolis-Hastings chain
    f    - function that implements the log unnormalized density of the desired equilibrium distribution
    init - initial state of the Markov Chain
    sigma - covariance matrix of the proposal distribution (default will be I)
    outfun - function that will transform the output of the chain

    returns a chain that uses the Metropolis-Hastings update
*)
let mhchain ?sigma ?(outfun=id) ~f init = 
  0

(* 

   Random-walk Metropolis algorithm. 
   f - normalized density
   init - initial state, a vector
   iterations - number of iterations to run the chain

   the proposal distribution is a multivariate normal with covariance matrix
   equal to the identity. it is symmetric and thus we use the metropolis update. 

 *)
let metropolis ?(rng=Distribution.stdrng) ~f init iterations = 
  let d = Matrix.vector_size init in    (* dimensionality *)
  let sigma = Matrix.identity d in
  let res = Matrix.create ~rows:iterations ~cols:d ~initval:0.0 in
  let rec iter i state = 
    Matrix.copy_vec_mat_row state res i;
    if i = iterations-1 then state
    else
      let newstate = Distribution.sample_mvnorm state sigma in
      let r = (f newstate) /. (f state) in
      if r >= 1.0 then iter (i+1) newstate
      else
        let u = Distribution.(rng.sample01 ()) in
        if u <= r then iter (i+1) newstate
        else iter (i+1) state in
  let _ = iter 0 init in
  res

let f1 state = 
  0.2

let () = 
  let initial = Matrix.vector_from_array [| 0.0; 0.0; 0.0 |] in
  let m = metropolis ~f:f1 initial 10 in
  Matrix.print m
