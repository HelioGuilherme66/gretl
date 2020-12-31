program define gretl_reorder
version 8.2
local matrix `1'
local fname `2'
tempname myfile
file open `myfile' using "`fname'", write text replace
local k = rowsof(`matrix')
local k1 = `k' - 1
file write `myfile' %8.0g (`k') %8.0g (`k') _n
file write `myfile' %15.0e (`matrix'[`k',`k'])
forvalues c=1/`k1' {
  file write `myfile' %15.0e (`matrix'[1,`c'])
}
file write `myfile' _n
forvalues r=1/`k1' {
  file write `myfile' %15.0e (`matrix'[`r',`k'])
  forvalues c=1/`k1' {
    file write `myfile' %15.0e (`matrix'[`r'+1,`c']) _n
  }
}
file close `myfile'
end
