method_labels <- c(
  siap = "SIAP",
  si = "SI",
  si_trend = "SI PrdMu",
  si_trend_cov = "SI PrdMu CrosSpec (Step 1)",
  sia = "SIA",
  sia_trend = "SIA LatentPrdMu (Step 2)",
  s1_si_trend = "Step 1 + SI LatentPrdMu",
  gp = "GP",
  marss = "MARSS",
  trmf = "TRMF",
  latc = "LATC"
)

method_colors <- c(
  siap = "#000000",      # black
  si = "blue",        
  si_trend = "lightgreen", 
  si_trend_cov = "cyan", 
  sia = "darkolivegreen",       
  sia_trend = "purple", 
  s1_si_trend = "darkslateblue", 
  gp = "red",        
  marss = "#F0E442",     
  trmf = "green",      
  latc = "darkgreen"  
)

method_shapes <- c(
  siap = 1,  # circle
  si = 2,    # triangle
  si_trend = 15, 
  si_trend_cov = 16, 
  sia = 17,       
  sia_trend = 18, 
  s1_si_trend = 5, 
  trmf = 0,  # square
  latc = 6,   # inverted triangle
  gp = 8,     # star
  marss = 4   # cross
)

realdata_colors <- c(csim = "purple", siap = "red",
            gp = "steelblue", marss = "gold2", tsis = "black", tsis_baseline = "darkgreen")

subset_methods <- c("siap", "si", "trmf", "latc", "gp", "marss")
siap_methods <- c("s1_si_trend", "si", "si_trend", "si_trend_cov", "sia", "sia_trend", "siap")
