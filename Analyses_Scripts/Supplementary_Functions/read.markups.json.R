# Function for importing landmark JSON files into R. Taken from https://github.com/SlicerMorph/SlicerMorphR
# Please cite the following when using this:
# Rolfe S., Pieper S., Porto A., Diamond K., Winchester J., Shan S., Kirveslahti H., Boyer D., Summers A., Maga A.M. 2021. SlicerMorph: An open and extensible platform to retrieve, visualize and analyse 3D morphology. Methods in Ecology and Evolution. 12:1816–1825.
# Rolfe S., Davis C., Maga A.M. 2021. Comparing semi-landmarking approaches for analyzing three-dimensional cranial morphology. American Journal of Physical Anthropology. 175:227–237.

read.markups.json = function(file=NULL){
  if (!require(jsonlite)) {
    print("installing jsonlite")
    install.packages('jsonlite')
    library(jsonlite)
  }
  dat = fromJSON(file, flatten=T)
  n=length(dat$markups$controlPoints[[1]]$position)
  labels = dat$markups$controlPoints[[1]]$label
  temp = array(dim = c(n, 3), dimnames=list(labels, c("X", "Y", "Z")))
  for (i in 1:n) temp[i,] = dat$markups$controlPoints[[1]]$position[[i]]
  return(temp)
}
