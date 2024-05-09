library(shiny)

library(survival)
library(survminer)
library(parallel)
library(promises)

ui <- fluidPage(
  "Hello, world!"
)
server <- function(input, output, session) {
}
shinyApp(ui, server)
