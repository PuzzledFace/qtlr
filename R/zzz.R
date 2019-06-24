#Hack to remove devtools'  warnings about no visible globals when using NSE.
Area <- NULL
Chain  <- NULL
Column  <- NULL
Dummy <- NULL
Events <- NULL
Index <- NULL
N <- NULL
ObservedResponse <- NULL
Parameter <- NULL
Sample <- NULL
Subjects <- NULL
Value <- NULL
a <- NULL
b <- NULL
x <- NULL
y <- NULL
Corps <- NULL
Deaths <- NULL
Event <- NULL
Guards <- NULL
Study <- NULL
SubjectID <- NULL
Year <- NULL
Group <- NULL
Censored <- NULL
Time <- NULL
TotalTime <- NULL

#' The Kalbfleisch and Prentice (1980) VA lung dataset.
#' 
#'  @format A tibble with 137 rows and 8 variables:
#'  \describe{
#'    \item{Treatment}{Treatment group}
#'    \item{CellType}{cell type}
#'    \item{SurvivalTime}{Survival time sice randomisation}
#'    \item{Status}{1=Event, 0=Censored}
#'    \item{KPS}{Karnofsky Performance Status}
#'    \item{Age}{Age in years at randomisation}
#'    \item{PriorTherapy}{Boolean.  received prior therapy?}
#'    }
#'  @source Kalbfleisch D and Prentice RL (1980), The Statistical Analysis of
#'  Failure Time Data. Wiley, New York 
"vaLung"

#' The Bortkiewicz cavalry dataset
#'
#' A dataset number of Prusiian cavalry soldiers kicked to death by their horses
#'
#' @format A tibble with 280 rows and 3 variables:
#' \describe{
#'   \item{Corps}{The corps}
#'   \item{Year}{The year}
#'   \item{Deaths}{The number of deaths in the given corps in the given year}
#'   ...
#' }
#' @source \url{https://archive.org/details/dasgesetzderklei00bortrich}
"cavalryDeaths"

#'Function to create the calavryDeaths tibble
makeCavalry <- function()
{
  cavalry <- utils::read.csv("~/Prussian-Horse-Kick-Data.csv")
  c <- tibble::as_tibble(cavalry) %>% 
         tidyr::gather(key=Corps, value=Deaths, Guards, tidyselect::starts_with("Corps")) %>% 
         dplyr::mutate(Corps=stringr::str_replace(Corps, "\\.", " "))
  cavalryDeaths <- c %>% dplyr::select(Corps, Year, Deaths)
  cavalryDeaths <- cavalryDeaths %>% 
                     dplyr::mutate(Corps=factor(Corps, labels=unique(cavalryDeaths$Corps), 
                                         levels=c("Guards", "Corps 1", "Corps 2", "Corps 3", "Corps 4", "Corps 5", "Corps 6", "Corps 7", "Corps 8", "Corps 9", "Corps 10", "Corps 11", "Corps 14", "Corps 15")))
  usethis::use_data(cavalryDeaths, overwrite=TRUE)
}

#Hack to to a quick and dirty reinstall during development
# reinstall <- function()
# {
#   utils::remove.packages("qtlr")
#   utils::install.packages(pkgs=devtools::build(), repos=NULL)
#   library("qtlr")
# }


