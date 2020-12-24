taxo_selected <- function(kingdom, supergroup, division, class){
  
# Small function to return the taxo level and taxon name 
# If the lowest level is "All" then it is the one selected then goes up one level

      if(class != "All") {
        taxo_level = "class"
        taxo_name = class
      } else {
        if(division != "All") {
        taxo_level = "division"
        taxo_name = division
        } else {
          if(supergroup != "All") {
          taxo_level = "supergroup"
          taxo_name = supergroup
          } else {
            taxo_level = "kingdom"
            taxo_name = kingdom          
          }
        }  
      }  
      
return( list(level = taxo_level, name = taxo_name))      
}

  print("taxonomy.R done")