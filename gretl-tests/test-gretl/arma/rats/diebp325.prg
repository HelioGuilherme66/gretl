*
* Sample random walks from pp 325-6.
*
all 300
set(first=1.0) y = y{1}+%ran(1.0)
set d = y-y{1}
*
* Because these are "random" numbers, they won't match the graph
* in the text, and, in fact, won't even match up if you run the
* program again. If you want to get the same values each time,
* you can use the SEED instruction at the top of the program.
* Just add SEED followed by a big integer (like SEED 12354). This
* works because computers actually generate what are known as
* pseudo-random numbers. There is actually a formula which generates
* them, but it has a very long cycle and produces numbers which pass
* just about any test for randomness.
*
graph(key=upleft,header='Figure 12.1 Random Walk: Level and Change') 2
# y
# d
*
* Random walk with drift
*
set(first=1.0) y = y{1}+.3+%ran(1.0)
set d = y-y{1}
graph(key=upleft,header='Figure 12.2 Random Walk with Drift: Level and Change') 2
# y
# d

