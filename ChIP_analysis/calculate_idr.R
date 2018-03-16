################################################################################
# Wrapper around idr package to calculate IDR for chip data
#
# Written by Michael Wolfe
#
# Copyright (c) 2018 Michael Wolfe University of Michigan. All rights reserved.
#
#
#Developed by: Michael Wolfe
#University of Michigan
#http://freddolino-lab.med.umich.edu/
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal with
#the Software without restriction, including without limitation the rights to
#use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
#the Software, and to permit persons to whom the Software is furnished to do so,
#subject to the following conditions:
#
#Redistributions of source code must retain the above copyright notice, this
#list of conditions and the following disclaimers.  Redistributions in binary
#form must reproduce the above copyright notice, this list of conditions and the
#following disclaimers in the documentation and/or other materials provided with
#the distribution.  Neither the names of Michael Wolfe, University of Michigan,
#nor the names of its contributors may be used to endorse or promote products
#derived from this Software without specific prior written permission.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
#FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE CONTRIBUTORS
#OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
#WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH THE SOFTWARE.
################################################################################
library(RcppCNPy)
library(ggplot2)
library(idr)
set.seed(1234) 
convert_to_robustz <- function(input_vec){
    finite <- is.finite(input_vec)
    this_median <- median(input_vec[finite])
    mad <- median(abs(input_vec[finite] - this_median))
    output <- (input_vec - this_median)/(mad*1.4826)
    output
}
    
args <- commandArgs(trailingOnly=TRUE)
act1 <- npyLoad(args[1])
act1 <- convert_to_robustz(act1)
act2 <- npyLoad(args[2])
act2 <- convert_to_robustz(act2)
actual <- cbind(act1, act2)
predicted_num <- as.numeric(args[3])
out_pre <- args[4]

mu <- 0.0
std <- 1.4826
rho <- 0.1
weight <- (predicted_num*250/10)/length(actual[,1])

this_data <- actual
this_data[is.na(this_data)] <- 0

print("Estimating IDR")
out.idr <- est.IDR(this_data, mu, std, rho, weight)
print("Done estimating IDR")
print(out.idr$para)
print(out.idr$loglik)
print("writing IDR")
png(paste0(out_pre, "_idr_llk.png"))
plot(out.idr$loglik.trace)
dev.off()
npySave(paste0(out_pre,"_idr.npy"), out.idr$IDR)
