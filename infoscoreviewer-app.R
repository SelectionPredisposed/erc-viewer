library(shiny)
library(dplyr)
library(pool)
library(RPostgreSQL)


u = Sys.getenv("POSTGRESQL_USER")
p = Sys.getenv("POSTGRESQL_PASSWORD")
print(p)
pool <- dbPool(
  drv = RPostgreSQL::PostgreSQL(max.con=20),
  dbname = "mobagenetics",
  host = "192.168.1.161",
  user = u,
  password = p
)
onStop(function() {
  poolClose(pool)
})

ui <- fluidPage(
  textInput("rsidin", "Enter your rsID:", "rs12"),
  DT::dataTableOutput("tbl")
)

server <- function(input, output, session) {
  # output$tbl <- DT::renderDataTable({
  #   pool %>% tbl("info_scores") %>% filter(rsid == input$rsid) %>% collect()
  # })

  output$tbl <- DT::renderDataTable({
    sql <- "SELECT * FROM info_scores WHERE rsid = ?rsidin;"
    #sql <- "SELECT * FROM info_scores WHERE rsid ~ ?rsidin::VARCHAR LIMIT 20;"
    #sql <- "SELECT * FROM info_scores WHERE ?rsidin::VARCHAR LIKE (rsid) OR ?rsidin::INTEGER LIKE (pos) LIMIT 20;"
    query <- sqlInterpolate(pool, sql, rsidin = input$rsidin)
    dbGetQuery(pool, query) %>% select(rsid, chr, pos, ref, alt, typed, info, ref_panel_af, info, cohort)
  })
}

shinyApp(ui, server)
