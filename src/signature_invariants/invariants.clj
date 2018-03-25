(ns
  signature-invariants.invariants
  "Generate the invariants of Diehl/Reizenstein - Invariants of multidimensional time series based on their iterated-integral signature, 2018."
  (:require [signature-invariants.young-tableau :as yt]
            [signature-invariants.signature :as signature]
            [signature-invariants.permutations :as permutations]
            [hopf-algebra.shuffle-algebra :as sa]
            [hopf-algebra.linear-combination :as lc]
            [clojure.math.combinatorics]
            [clatrix.core :as clatrix]
            ))

;;;;;;;;;;;;;;;;;;
; HELPER FUNCTIONS
;;;;;;;;;;;;;;;;;;

(defn index-generator
  "Returns a list of all sequences of length `length` in the numbers 0,..,`n-letters`-1."
  [length n-letters]
  (apply clojure.math.combinatorics/cartesian-product (repeat length (range n-letters))))

(defn identiy-matrix [dim]
  (vec (repeat dim (vec (repeat dim 0)))))

(defn rank
  "The linear rank of the collection of vectors `coll`."
  [coll]
  (if (empty? coll)
    0
    (clatrix/rank (clatrix/matrix coll))))

;;;;;;;;;;;;;;;
; GL INVARIANTS
;;;;;;;;;;;;;;;

; First, a brute force calculation, as a sanity check.

(defn- g-brute-force [i perm dim]
  (let [weight (/ (count perm) dim)]
    (reduce
      *
      (for [z (range weight)]
        (permutations/det
          (loop [m (identiy-matrix dim)
                 j 0]
            (if (= j dim)
              m
              (recur
                (assoc-in m [ (nth i (nth perm (+ (* z dim) j))) j] 1)
                (inc j)))))))))

(defn- gl-invariants-brute-force-for-permutation [dim weight perm]
  (loop [result {}
         indices (index-generator (* dim weight) dim)]
    (if (empty? indices)
      result
      (recur
        (lc/lc-add
          result
          (lc/lc-multiply
            (g-brute-force (first indices) perm dim) ; XXX most of these are zero ..
            { (sa/->ShuffleWord (first indices) ) 1 }))
        (rest indices)))))

(defn gl-invariants-brute-force
  "Returns a generating set for GL(R^dim) invariants of weight `weight`"
  [dim weight]
  (for [perm (clojure.math.combinatorics/permutations (range (* dim weight)))]
    (gl-invariants-brute-force-for-permutation dim weight perm)))


; Now, the efficient way using Young tableaux.

(defn permute
  "i_{\\sigma(1) .. i_{\\sigma(n)}.
  Here, `sigma` is either a map of the form {0 0, 1 2, 2 1},
  or a list in one row notation, in this example [0 2 1]."
  [sigma i]
  (if (= hopf_algebra.shuffle_algebra.ShuffleWord (type i)) ; XXX ugly case destinction
    (sa/->ShuffleWord (permute sigma (:content i)))
    (let [sigma'
          (if (map? sigma)
            sigma
            (into {} (map vector (range (count sigma)) sigma)))]
    (for [r (range (count i))]
      (nth i (get sigma' r r))))))

(defn- overall-sign [ taus-with-signs ]
  (reduce * (map second taus-with-signs)))

(defn- index-for-taus [ dim weight taus-with-signs ]
  (reduce
    (fn [i t] (permute t i))
    (apply concat (repeat weight (range dim)))
    (map first taus-with-signs)))

(defn- taus
  "Returns list of tuples of permutations and their sign.
  For example, for `dim`=3, `weight`=2 in each tuple are two permutations
  the first one of [0 1 2], the second on of [3 4 5]."
  [dim weight]
  (apply
    clojure.math.combinatorics/cartesian-product
    (map
      (fn [w] 
        (map
          (fn [ [perm sign] ]
            [ (into {} (map vector (range (* w dim) (* (inc w) dim))
                            (map (fn [x] (+ (* dim w) x)) perm))),
               sign] )
          (permutations/signed-permutations dim)))
      (range weight))))

(defn- gl-invariant-for-identity [dim weight]
  (loop [result {}
         taus-with-signs-s (taus dim weight)]
    (if (empty? taus-with-signs-s)
      result
      (recur
        (lc/lc-add
          result
          (lc/lc-multiply
            (overall-sign (first taus-with-signs-s))
            { (sa/->ShuffleWord (index-for-taus dim weight (first taus-with-signs-s))) 1 })) ; XXX sigma vs sigma^-1
        (rest taus-with-signs-s)))))

(defn to-mathematica
  "Output a lc of words to an array that can be used in Mathematica."
  [lc]
  (println "{")
  (doseq [ [k v] lc ]
    (println "{ {" (apply str (interpose "," (map inc (:content k)))) "}, " v"},"))
  (println "}"))

(defn gl-invariants
  "Returns a basis for GL(R^dim) invariants of weight `weight`"
  [dim weight]
  (let [for-identity (gl-invariant-for-identity dim weight)
        sigmas (map yt/tableau->permutation
                    (yt/standard-tableaux (repeat dim weight)))]
    (for [sigma sigmas]
      (do
      (lc/lc-apply-linear-function
        (fn [i] { (permute (yt/invert sigma) i) 1}) ; XXX sigma vs sigma^-1
        for-identity)))))

;;;;;;;;;;;;;;;
; SO INVARIANTS
;;;;;;;;;;;;;;;

(defn- even-number-of-1s [order]
  ; XXX slow hack
  (filter
    (fn [i] (even? (count (filter #( = % 1 ) i))))
    (index-generator order 2))
  )
(defn- odd-number-of-1s [order]
  ; XXX slow hack
  (filter
    (fn [i] (odd? (count (filter #( = % 1 ) i))))
    (index-generator order 2))
  )

(defn- same-number-of-0s-1s [order]
  ; XXX slow hack
  (filter
    (fn [i] (=
             (count (filter #( = % 0 ) i))
             (count (filter #( = % 1 ) i))))
    (index-generator order 2))
  )

(defn- exp [x n]
  (reduce * (repeat n x)))

(defn- real [i]
  (let [order (count i)]
    (reduce
      lc/lc-add
    (for [j (even-number-of-1s order)]
      (lc/lc-multiply
        (*
          (exp -1 (/ (count
                    (filter
                      (fn [ [a b] ] (= 1 a))
                      (map vector j i))) 2))
          (exp -1 (count
                    (filter
                      (fn [ [a b] ] (and (= 1 a) (= 1 b)))
                      (map vector j i)))))
        {(sa/->ShuffleWord j) 1})))))

(defn- im [i]
  (let [order (count i)]
    (reduce
      lc/lc-add
    (for [j (odd-number-of-1s order)]
      (lc/lc-multiply
        (*
          (exp -1 (/ (dec (count
                    (filter
                      (fn [ [a b] ] (= 1 a))
                      (map vector j i)))) 2))
          (exp -1 (count
                    (filter
                      (fn [ [a b] ] (and (= 1 a) (= 1 b)))
                      (map vector j i)))))
        { (sa/->ShuffleWord j) 1})))))

(defn so2-invariants [order]
  (assert even? order)
  (flatten
    (for [i (same-number-of-0s-1s order) ; XXX very slow
          :when (= (first i) 0)] 
    [(real i) (im i)])))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; jeremy's algorithm https://github.com/bottler/iisignature/blob/master/src/rotationalInvariants.hpp

;  //Returns all the combinations of n 0s and n 1s which begin with 0
;  vector<vector<unsigned char>> possibleHalves(int n) {
;    vector<unsigned char> base(n + n, 0);
;    std::fill(base.begin() + n, base.end(), (unsigned char) 1);
;    vector<vector<unsigned char>> out{ base };
;    while (std::next_permutation(base.begin(), base.end())) {
;      if (base[0] == 1)
;        break;
;      out.push_back(base);
;    }
;    return out;
;  }

;(defn possible-halves
;  "Returns all the combinations of order/2 0s and order/2 1s which begin with 0"
;  [order]
;  (let [n (/ order 2)
;        start-vector (concat (repeat (dec n) 0) (repeat n 1))]
;    (for [perm (clojure.math.combinatorics/permutations start-vector)]
;          (cons 0 perm))))
;
;(defn multiply-out-terms
;  "expression being {1,0,0,1} means (x+iy)(x-iy)(x-iy)(x+iy)
;   this returns the real and imaginary parts"
;  [expression]
;  (loop [real '( (0, 1) )
;         imag '( )
;         remaining-expression expression
;         ]
;    (if (empty? expression)
;      [real imag]
;      ; TODO continue here
;      (let []
;        (recur 
;          )))))

;;  //an Idx represents a word in two variables as a binary number
;;  //- and therefore as an offset into a level in the signature.
;;  // (If you set start1, it will be preceeded by a 1, so e.g. "112" 
;;  // and "12" will be represented as 1001b and 101b respectively
;;  // instead of both being 1, which may be easier to keep track of.
;;  // The rest of the code would have to change though.)
;;  using Idx = std::uint64_t;
;;  using Invariant = vector<pair<Idx, double>>;

;; //expression being {1,0,0,1} means (x+iy)(x-iy)(x-iy)(x+iy)
;;   //this returns the real and imaginary parts 
;;   pair<Invariant, Invariant> multiplyOutTerm(vector<unsigned char>& expression)
;;   {
;;     bool start1 = false;
;;     Invariant real, imag;
;;     real.reserve(((size_t)1u) << expression.size());
;;     imag.reserve(((size_t)1u) << expression.size());
;;     real.emplace_back(start1 ? 1 : 0, 1);
;;     for (unsigned char c : expression)
;;     {
;; #if 0 
;;       //Some unfortunate linux setups have new GCC but a std library so old that this
;;       //version of insert is a void function
;;       auto realHalf = real.insert(real.end(), imag.begin(), imag.end());
;;       auto imagHalf = imag.insert(imag.end(), real.begin(), realHalf);
;; #else
;;       auto realSize = real.size(), imagSize = imag.size();
;;       real.insert(real.end(), imag.begin(), imag.end());
;;       auto realHalf = real.begin() + realSize;
;;       imag.insert(imag.end(), real.begin(), realHalf);
;;       auto imagHalf = imag.begin() + imagSize;
;; #endif
;;       std::for_each(real.begin(), realHalf, [](pair<Idx, double>& p) {
;;         p.first *= 2;
;;       });
;;       std::for_each(imag.begin(), imagHalf, [](pair<Idx, double>& p) {
;;         p.first *= 2;
;;       });
;;       std::for_each(realHalf, real.end(), [c](pair<Idx, double>& p) {
;;         p.first = p.first * 2 + 1;
;;         if (!c)
;;           p.second *= -1;
;;       });
;;       std::for_each(imagHalf, imag.end(), [c](pair<Idx, double>& p) {
;;         p.first = p.first * 2 + 1;
;;         if (c)
;;           p.second *= -1;
;;       });
;;     }
;;     return { std::move(real), std::move(imag) };
;;   }

;(defn so2-invariants-jeremy [order]
;  (assert even? order)
;
;  )

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defn- equal-size-partitions
  "Returns all partitions of 0,..,`n`-1
   into set of size `size`."
  [n size]
  (filter
    (fn [x] (every? #( = size (count %)) x)) ; XXX this is expensive
    (clojure.math.combinatorics/partitions (range n))))

(defn- multiply-specific-positions
  " order = 3
    positions-monomials = [ [ [0 3] [7 1] ], [ [2] [3] ] ]
   
    => [7 3 1]"
  [order positions-monomials]
  (loop [result (vec (repeat order 0))
         positions-monomials positions-monomials]
    (if (empty? positions-monomials)
      result
      (recur
        (let [ [p m] (first positions-monomials) ]
          (apply assoc result (interleave p m)))
        (rest positions-monomials)))))

(defn- lc-multiply-specific-positions
  [order positions lcs]
  (loop [result {}
         remaining (apply clojure.math.combinatorics/cartesian-product lcs)]
    (if (empty? remaining)
      result
      (recur
        (let [k-v (first remaining)]
        (assoc
          result
          (multiply-specific-positions order (map vector positions (map key k-v)))
          (reduce * (map val k-v))))
        (rest remaining)))))

(defn- so-invariants-without-det [dim order]
  (if (odd? order)
    []
    (let [inner (into {} (map (fn [x] [ [x x] 1 ]) (range dim)))]
      (for [part (equal-size-partitions order 2)]
        (lc-multiply-specific-positions order
                                        part
                                        (repeat (/ order 2) inner))))))

(defn- subsets [coll size]
  (clojure.math.combinatorics/combinations coll size))

(defn- det-helper [dim]
  (into {}
    (for [ [perm sign] (permutations/signed-permutations dim)]
      [ perm sign ])))

(defn- so-invariants-with-one-det [dim order]
  (if (or (< order dim) ; no determinant fits
          (odd? (- order dim))) ; using a determinant, the others cannot be used for inner products
    []
    (let [indices (set (range order)),
          det (det-helper dim)
          tmp (so-invariants-without-det dim (- order dim))]
      (for [subset (subsets indices dim),
            tmpp tmp,
            :let [subset-complement (clojure.set/difference indices subset)] ]
        (if (empty? subset-complement)
          det
          (lc-multiply-specific-positions
            order
            [subset subset-complement]
            [det tmpp]))))))

(defn so-invariants
  "Returns a generating set for SO(R^dim) invariants of order `order`"
  [dim order]
  (map
    (fn [lc] (lc/lc-apply-linear-function (fn [x] {(sa/->ShuffleWord x) 1}) lc)) ; wrap in ShuffleWord
    (concat
      (so-invariants-without-det dim order)
      (so-invariants-with-one-det dim order))))


;;;;;;;;;;;;;;;;;;;;;;;;
; PERMUTATION INVARIANTS
;;;;;;;;;;;;;;;;;;;;;;;;

(defn- nabla [ell]
  ; XXX this is a very expensive brute force implementation
  (into #{}
    (let [letters (into #{} ell)]
      (for [s letters]
        (map first (filter (fn [[i x]] (= x s)) (map-indexed vector ell)))))))

(defn permutation-invariants
  "Returns a basis for the permuation invariants of order `order` and dimension `dim`."
 [dim order]
 (let [parts (clojure.math.combinatorics/partitions (range order))]
   (for [p (map (fn [x] (into #{} x)) parts)]
     (loop [result {}
            remaining-indices (index-generator order dim)]
       (if (empty? remaining-indices)
         result
         (recur
           (if (= p (nabla (first remaining-indices)))
             (lc/lc-add
               {(sa/->ShuffleWord (first remaining-indices)) 1}
               result)
             result)
             (rest remaining-indices)))))))

;;;;;;;;;;;;;;;
; GENERAL STUFF
;;;;;;;;;;;;;;;

(defn- zeroes [n]
   (vec (take n (repeat 0))))

(defn count-independents
  "`coll` is a collection of lc's containing ShuffleWord/ConcatWord
    in dimension `dim` of order `order`.
    The dimension of the subspace spanned is returned."
  [coll dim order]
  (if (empty? coll)
    0
    (let [matrix (for [x coll] (signature/as-vector x dim order))]
      (rank matrix))))

(defn just-independents
  "`coll` is a collection of lc's containing ShuffleWord/ConcatWord
    in dimension `dim` of order `order`.
    A subset of maximal linearly independents is returned."
  [coll dim order]
   ; XXX seems way too complicated
  (loop [
          rank-so-far 0
          so-far-matrix []
          so-far-raw []
          remaining coll ]
    (if (empty? remaining)
      so-far-raw
      (let [next-vector (signature/as-vector (first remaining) dim order)
            new-rank (rank (conj so-far-matrix next-vector))
            add-it?  (or (> new-rank rank-so-far))]
      (recur
        new-rank
        (if add-it? (conj so-far-matrix next-vector) so-far-matrix)
        (if add-it? (conj so-far-raw (first remaining)) so-far-raw)
        (rest remaining))))))
