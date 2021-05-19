# --- Supergroups
supergroups <- reactive({
  req(input$kingdom)
  filter(pr2_taxo, kingdom == input$kingdom) %>% 
    pull(supergroup) %>%  
    unique()  
})

observeEvent(supergroups(), {
  updateSelectInput(session, "supergroup", choices = c("All", supergroups()))
})  

# --- Divisions
divisions <- reactive({
  req(input$supergroup)
  filter(pr2_taxo, supergroup == input$supergroup) %>% 
    pull(division) %>%  
    unique()  
})

observeEvent(divisions(), {
  updateSelectInput(session, "division", choices = c("All", divisions()))
})

# --- Classes    
classes <- reactive({
  req(input$division)
  filter(pr2_taxo, division == input$division) %>% 
    pull(class) %>%  
    unique()
})

observeEvent(classes(), {
  updateSelectInput(session, "class", choices = c("All", classes()))
})