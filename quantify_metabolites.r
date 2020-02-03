
library(argparse)
#install.packages('speaq')
#library(speaq)


parser <- ArgumentParser()

parser$add_argument('-dir', type='character', required=TRUE)
parser$add_argument('-csv', type='character', required=TRUE)
parser$add_argument('-out', type='character', required=TRUE)
parser$add_argument('-ncores', type='integer', default=1)
parser$add_argument('-library', type='character', default='default')


args = parser$parse_args()

library(ASICS)
library(ASICSdata)

if (args$library != 'default') {
  Flib = scan(args$library)
  used_library <- pure_library[c(Flib)]
} else {
  used_library <- pure_library
}


spectra_data <- importSpectra(
  name.dir = args$dir,
  name.file = args$csv,
  type.import = 'csv'
)

spectra_align <- alignSpectra(spectra_data)
spectra_obj <- createSpectra(spectra_align)

ASICS_results <- ASICS(spectra_obj, ncores = args$ncores, pure.library = used_library)

quant <- getQuantification(ASICS_results)

write.csv(quant, file = args$out)