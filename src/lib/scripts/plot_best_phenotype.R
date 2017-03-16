
#***************************************************************************
# Copyright (C) 2014-2017 Charles Rocabert, Carole Knibbe, Guillaume Beslon
# Web: https://github.com/charlesrocabert/Evo2Sim
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#***************************************************************************

### generate the list of prime numbers ###
# array indexes correspond to numbers.
# if value is 1, the number is prime number.
# (else 0)
##########################################
build_prime_numbers_list <- function( maximum )
{
  prime_numbers = rep(1, maximum)
  prime_numbers[1] = 0;
  for (i in seq(2, maximum))
  {
    multiple = 2*i
    while (multiple <= maximum)
    {
      prime_numbers[multiple] = 0
      multiple = multiple+i
    }
  }
  return(prime_numbers)
}

#-------------------------------------#
# 1) Read command line arguments      #
#-------------------------------------#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0)
{
  stop("You must give the experiment path")
}
experiment_path = args[1]
setwd(experiment_path)

silence_warnings = F
if (length(args) == 2)
{
  if (args[2] == "-no-warnings")
  {
    silence_warnings = T
  }
}

#-------------------------------------#
# 2) Hide error messages and warnings #
#-------------------------------------#
if (silence_warnings == T)
{
  options(show.error.messages = FALSE, warn=-1)
}

#-------------------------------------#
# 3) Open data                        #
#-------------------------------------#
d = read.table("./statistics/best_metabolic_amounts.txt", sep=" ")
current = d[d[,2]>0,]
inherit = d[d[,3]>0,]

#-------------------------------------#
# 4) Plot the state vector            #
#-------------------------------------#
prime_numbers = build_prime_numbers_list(5000)
colors = c()
for (i in seq(1, length(current[,1])))
{
  if (prime_numbers[current[i,1]] == 0)
  {
    colors = c(colors, "darkolivegreen3")
  }
  if (prime_numbers[current[i,1]] == 1)
  {
    colors = c(colors, "tomato")
  }
}

png("./figures/best_metabolic_vector.png", width=1200, height=600, pointsize=20)
barplot(current[,2], names.arg=current[,1], col=colors, main="Best metabolic state vector", cex.names=0.5, xlab="Metabolite tag", ylab="Concentration")
barplot(inherit[,3], border="black", col=rgb(0,0,0,0.25), add=T)
legend("topright", legend=c("Essential", "Non essential", "Inherited"), col=c("tomato", "darkolivegreen3", "black"), lwd=c(2,2))
dev.off()

