*
* Sample autocorrelation graphs from pp 121-123
*
all 36
*
* It's a good idea to graph autocorrelation functions with the
* options min=-1.0,max=1.0, so they will always have the same
* vertical scale. A bar graph is the most common choice for the
* style. The number=0 option makes the axis numbering start with
* 0.
*
set graddamp = .9**t
graph(header='Figure 6.2 Autocorrelation Function, Gradual One-Sided Damping',number=0,style=bar,max=1.0,min=-1.0)
# graddamp
*
set nodamp = %if(t==1,1.0,.95)
graph(header='Figure 6.3 Autocorrelation Function, No Damping',number=0,style=bar,max=1.0,min=-1.0)
# nodamp
*
set osdamp = .9**(t-1)*(cos((t-1)))
graph(header='Figure 6.4 Autocorrelation Function, Gradual Damped Oscillation',number=0,style=bar,max=1.0,min=-1.0)
# osdamp
*
set cutoff = %if(t<=14,.97**(t-1),0.0)
graph(header='Figure 6.5 Autocorrelation Function, Sharp Cutoff',number=0,style=bar,max=1.0,min=-1.0)
# cutoff

