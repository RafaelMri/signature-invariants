# signature-invariants

A collection of functions that calculate the invariants
of Diehl/Reizenstein - Invariants of multidimensional signals based on their signature, 2018.

For example usage, see the worksheets folder.


## Running the code

Installation is simple if you have [leiningen](http://leiningen.org); this tool will arrange
to retrieve everything else you need. On Mac OS, for example,

~~~ sh
$ brew install leiningen
~~~

ought to get you started. Clone the repo. In the repo folder, `lein repl` will start a REPL,
where you can run the examples from above.

To run the test suite:

~~~ sh
$ lein test
~~~

To run the worksheets:

~~~ sh
$ lein gorilla
~~~
Then open the output URL in your browser and press <CTRL-G> <CTRL-G>, and choose worksheets/worksheet-invariants.clj.
You can then run the sheet with <CTRL> <SHIFT> <ENTER>.

Copyright © 2018 Joscha Diehl

Distributed under the [Eclipse Public License](https://opensource.org/licenses/eclipse-1.0.php), the same as Clojure.
