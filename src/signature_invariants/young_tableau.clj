(ns signature-invariants.young-tableau
  "Methods for generating and manipulating Young tableaux.")

(defn first-standard-tableau [shape]
  "the first standard tableau in the order used here;
   this means filling _columns_ first"
  (let [shape (vec shape)
        size (apply + (flatten shape))]
    (loop [initial-tableau (vec (for [s shape] (vec (range s))))
                     row 0
                     col 0
                     i 0]
                (if (= i size)
                  initial-tableau
                  (let [next-row? (and (< (inc row) (count shape)) (< col (shape (inc row))))]
                    (recur
                      (assoc-in initial-tableau [row col] i)
                      (if next-row?
                        (inc row)
                        0)
                      (if next-row?
                        col
                        (inc col))
                      (inc i)))))))

(defn positions-for-f [f coll]
  (keep-indexed (fn [idx val] (when (f val) idx)) coll))

(defn positions [x coll]
  (positions-for-f (fn [z] (= z x)) coll))

(defn standard-tableau->vector [tableau]
  " [0 1 3]
    [2]         -> [0 0 1 0 2]
    [4]  "
  {:pre [(= (into #{} (flatten tableau))
            (into #{} (range (apply + (map count tableau)))))]}
  (let [size (apply + (map count tableau))]
    (loop [i 0
           v (vec (repeat size 0))]
      (if (= i size)
        v
        (recur
          (inc i)
          (assoc v i
                 (first (positions-for-f (fn [row] (some #(= i %) row)) tableau))))))))

(defn vector->standard-tableau [v]
  (vec (for [row (range (inc (apply max v)))] (vec (positions row v)))))

(defn last-standard-tableau [shape]
  "the first standard tableau in the order used here;
   this means filling _rows_ first"
  (vector->standard-tableau
    (apply concat
           (map-indexed (fn [ i s ] (vec (repeat s i))) shape))))

(defn- find-j [v]
  (loop [j 1
         ell (-> (vec (repeat (count v) 0))
                 (update (nth v 0) inc)
                 (update (nth v 1) inc))]
    (if (< (nth v j) (nth v (dec j)))
      [j ell]
      (recur
        (inc j)
        (update ell (nth v (inc j)) inc)))))

(defn rows->cols [rows]
  (loop [rows rows
         cols [] ]
    (if (= 0 (first rows))
      cols
      (recur
        (map (fn [x] (if (> x 0) (dec x) x)) rows)
        (conj cols (count (filter #( > % 0) rows)))))))

(defn- change-j [v j ell]
  (let [j-row (nth v j),
        row-below (nth ell (inc j-row))
        target (last (positions row-below ell))]
    [ (assoc v j target) (update ell target dec) ]))

(defn- fill-T-j [v j ell]
  (vec
    (concat
      (apply concat (map (comp vec range) (rows->cols ell)))
      (subvec v j))))

;;
;- the algorithm from StandardTableaux_shape.__iter__ in SAGE
;;
(defn- standard-tableaux-ntv [current-vector]
  (let [ ;- find first j s.t. j is _not_ in the lowest corner of T_j := 1,..,j; ell = shape of T_j
         [j ell] (find-j current-vector)
         ;- put j in the right corner of the lowest row of T_j that has same size as the row below j
         [next-vector ell'] (change-j current-vector j ell)
         ;- fill in the rest of T_j with !,..,j-1; filling the _columns_ first
         next-vector' (fill-T-j next-vector j ell')]
    next-vector'
  ))

(defn standard-tableaux [ shape ]
  (let [size (apply + shape)
        first-one (first-standard-tableau shape)
        last-one (last-standard-tableau shape)
        ]
    (loop [current-tableau-vector (standard-tableau->vector first-one)
           result [ first-one ] ]
      (if (= last-one (last result))
        result
        (let [next-tableau-vector (standard-tableaux-ntv current-tableau-vector)]
          (recur
            next-tableau-vector
            (conj result (vector->standard-tableau next-tableau-vector))))))))

(defn tableau-rows->tableau-cols [rows]
  (loop [rows rows
         cols [] ]
    (if (empty? rows)
      cols
      (recur
        (remove empty? (map rest rows))
        (conj cols (map first rows))))))

(defn is-standard-tableau? [t]
  (and
    (= (sort (flatten t)) (range (count (flatten t))))
    (every?  (fn [row] (= (sort row) row)) t)
    (every?  (fn [col] (= (sort col) col)) (tableau-rows->tableau-cols t))))

(defn invert [tau]
  (into {} (map (fn [ [k v] ] [v k]) tau))
  )

(defn tableau->permutation [t]
  "the permutation \\sigma such that \\sigma initial-tableau = t"
  ;(println "t="t)
  (let [initital-tableau (first-standard-tableau (mapv count t))]
    (into {} (map vector (flatten initital-tableau) (flatten t)))))

