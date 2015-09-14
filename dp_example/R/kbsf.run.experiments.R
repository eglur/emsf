source("kbsf.experiments.R")


puddle.experiment.kbrl.equi <- function(num.archs = seq(50, 250, by=50)*4) {
   pe <- load.sars("./kbsf/puddle/equi/puddle_data")
   puddle.kbrl(make.lsars(pe), filename = "./kbsf/puddle/equi/puddle_kbrl_")
     }


puddle.experiment.kbrl.uni <- function() {
   pu <- load.sars("./kbsf/puddle/uni/puddle_data")
   puddle.kbrl(make.lsars(pu), filename = "./kbsf/puddle/uni/puddle_kbrl_")
     }


puddle.experiment.kbrl.exp <- function(num.archs = seq(50, 250, by=50)*4) {
   pe <- load.sars("./kbsf/puddle/exp/puddle_data")
   puddle.kbrl(make.lsars(pe), filename = "./kbsf/puddle/exp/puddle_kbrl_")
   }


puddle.experiment.lspi.kbsf.equi<- function(num.archs = seq(50, 250, by=50)*4){
   pe <- load.sars("./kbsf/puddle/equi/puddle_data")
   puddle.lspi.kbsf(pe, filename.lspi ="./kbsf/puddle/equi/puddle_lspi_",
                        filename.kbsf="./kbsf/puddle/equi/puddle_kbsf_")
     }


puddle.experiment.lspi.kbsf.uni <- function() {
   pu <- load.sars("./kbsf/puddle/uni/puddle_data")
   puddle.lspi.kbsf(pu, filename.lspi ="./kbsf/puddle/uni/puddle_lspi_",
                        filename.kbsf="./kbsf/puddle/uni/puddle_kbsf_")
     }


puddle.experiment.lspi.kbsf.exp <- function(num.archs = seq(50, 250, by=50)*4) {
   pe <- load.sars("./kbsf/puddle/exp/puddle_data")
   puddle.lspi.kbsf(pe, filename.lspi ="./kbsf/puddle/exp/puddle_lspi_",
                        filename.kbsf="./kbsf/puddle/exp/puddle_kbsf_")
   }
