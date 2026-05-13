#Name: Sumed Seeyakmani Kuson


library(glmnet)
##CHILD Score, CHILD Stage, MELD, Age, BL_DACL, Albumin, CRP
library(survival)
library(openxlsx)
library(dplyr)
library(tidyr)
library(riskRegression)
library(pec)
library(randomForestSRC)
df <- vienna2_IRSQR_1005_surv
print(colnames(df), max = 1e6)   
cat(paste(sprintf('"%s"', colnames(df)), collapse = ", "), "\n")
library(rpart.plot)

str(df)

## =========================
## 1) MISSINGNESS SUMMARY
## =========================

# Counts how many rows have at least ONE missing value (NA) anywhere in the row.
rows_with_na <- sum(!complete.cases(df))
cat("Rows with at least one NA:", rows_with_na, "\n")

# Counts how many NAs are in each column.
na_summary <- sapply(df, function(x) sum(is.na(x)))

# Put the NA counts into a neat table for viewing/exporting.
na_table <- data.frame(
  Column = names(na_summary),
  NAs = as.integer(na_summary)
)

# Adds % missing per column (good for reporting / deciding what to drop/impute).
na_table$Percent <- round(na_table$NAs / nrow(df) * 100, 2)

# Exports the NA summary table to Excel for documentation / sharing.
write.xlsx(
  na_table,
  file = "C:/Users/Sumed Seeyakmani/Downloads/NA_summary.xlsx",
  rowNames = FALSE
)

cat(" NA summary exported successfully to 'NA_summary.xlsx' in your working directory.\n")

## =========================
## 2) CLEANUP / PREP
## =========================

# Drops columns we don't want in the analysis (e.g., list-columns or unused blocks).
df <- df %>% select(-c(clinical_characteristics, lab_data, baprof, baprof_exploratory_variables, fu_events))

# Ensures survival time is strictly > 0 (Cox model cannot handle time <= 0).
df$tfs_months <- ifelse(df$tfs_months <= 0, 1e-6, df$tfs_months)

# Ensures event coding is integer (helpful for consistency / modeling functions).
df$tfs_event_012 <- as.integer(df$tfs_event_012)

# Declares childstage as categorical (factor) so Cox treats it as a group variable.
df$childstage <- as.factor(df$childstage)

# Ensures the baseline dACLD variable has the intended reference level ordering.
# (First level "cACLD" becomes the reference in regression outputs.)
df$bl_dacld_01 <- factor(df$bl_dacld_01, levels = c("cACLD", "dACLD"))


## =========================
## 3) CREATE TWO EVENT DEFINITIONS
## =========================

# Event 1: "death-as-event" (tfs_event_012 == 2 means event happened).
# Everything else (0/1) is censored in this definition.
df$event_death <- ifelse(df$tfs_event_012 == 2L, 1L, 0L)

# Event 2: "transplant-as-event" (tfs_event_012 == 1 means event happened).
# Everything else (0/2) is censored in this definition.
df$event_transplant <- ifelse(df$tfs_event_012 == 1L, 1L, 0L)


# Example on how Cox model is implimented
# ## =========================
# ## 4) COX MODEL TEMPLATE
# ## =========================
# 
# # A Cox model template with a placeholder for the event column.
# # This lets us reuse the same predictor set but swap the event definition.
# cox_form <- as.formula(
#   Surv(tfs_months, EVENT_PLACEHOLDER) ~ childstage + meld + age + bl_dacld_01 + alb + crp
# )
# 
# 
# ## =========================
# ## 5) HELPER FUNCTION: FIT + REPORT + PH TEST
# ## =========================
# 
# fit_report_ph <- function(data, event_col) {
#   
#   # Converts the formula to text so we can replace EVENT_PLACEHOLDER with the real event column name.
#   ftxt <- deparse(cox_form)
#   ftxt <- sub("EVENT_PLACEHOLDER", event_col, ftxt, fixed = TRUE)
#   fml  <- as.formula(ftxt)
#   
#   # Extracts all variable names needed by the model (time, event, predictors).
#   # Then keeps ONLY complete cases across those vars (Cox/coxph needs no NA rows).
#   vars_needed <- all.vars(fml)
#   dsub <- data[complete.cases(data[, vars_needed]), vars_needed]
#   
#   # Safety: ensure time > 0 inside the model dataset.
#   
#   dsub$tfs_months <- ifelse(dsub$tfs_months <= 0, 1e-06, dsub$tfs_months)
#   
#   # Simple logging: which event, and how many rows used after dropping NAs.
#   cat("\n============================================================\n")
#   cat("Fit Cox model with event column:", event_col, "\n")
#   cat("N rows used:", nrow(dsub), "\n")
#   
#   # Fits the Cox proportional hazards model (Efron is standard for handling ties).
#   fit <- coxph(fml, data = dsub, ties = "efron")
#   print(summary(fit))  # shows HRs, SEs, p-values, concordance, etc.
#   
#   # Tests the proportional hazards assumption using Schoenfeld residuals.
#   # If a covariate has small p-value, PH may be violated for that covariate.
#   cat("\n--- Proportional Hazards (Schoenfeld) test ---\n")
#   ph <- cox.zph(fit)
#   print(ph)
#   
#   # Returns the fitted model object invisibly (so we can store it without printing again).
#   invisible(fit)
# }

# 
# ## =========================
# ## 6) RUN COX MODELS FOR BOTH EVENT DEFINITIONS
# ## =========================
# 
# # Fit Cox with death-as-event coding.
# fit_death <- fit_report_ph(df, event_col = "event_death")
# 
# # Fit Cox with transplant-as-event coding.
# fit_tx <- fit_report_ph(df, event_col = "event_transplant")


# ## =========================
# ## 7) CONCORDANCE (DISCRIMINATION)
# ## =========================
# 
# # Concordance index (C-index): how well the model ranks risk.
# # 0.5 = random, 1.0 = perfect ranking (in survival sense).
# cat("\nConcordance (death-as-event):    ", summary(fit_death)$concordance[1], "\n")
# cat("Concordance (transplant-as-event):", summary(fit_tx)$concordance[1], "\n")
# 

## =========================
## 8) QUICK CROSS-TAB CHECK
## =========================

# Pulls out childstage and childscore, then tabulates their joint counts.
# Useful to see if stage is consistent with score categories, check coding issues, etc.
df1 <- df$childstage
df2 <- df$childscore
table(df1, df2)


##LASSO REGRESSION

# 1st do Lasso for clinical variables
#"hvpg_2gr_16"                  "bl_stiffness"                
#"icg-pdr"                      "pcg-r15"                    
# From the Clinical Characteristics, hvpg_2gr_16, blstiffness, icg and pcg wont be used in the regression.
# 22 Variables
clinical_vars <- c(
  "sex_code", "age", "hvpg_strata_3gr_16", "hvpg_strata_3gr_20",
   "bl_easl_5_multistate", "bl_dacld_01", "varbleed_01",
  "asc_01", "he_01", "hcc1_bl01", "etio_systematic", "etio_code",
   "varices_simple", "bl_asc_123", "bl_he_123",
  "childscore", "childstage", "meld", "cni_01", "bl_dm_012",
  "bl_bmi", "hvpg"
)


# Example on how to use glmnet for LASSO regularization
## =========================
## 9) LASSO COX (GLMNET) — CLINICAL VARS
## =========================


# df_clinical_filtered <- df[complete.cases(df[, clinical_vars]), ]
# 
# # Print how many rows remain and how many were removed
# cat("Rows remaining :", nrow(df_clinical_filtered),"\n")
# 
# 
# set.seed(123)  # for reproducibility
# 
# # Define event and predictors 
# y <- Surv(df$tfs_months, df$event_death)
# vars <- clinical_vars
# 
# # Keep only rows with no NA in these variables and survival columns
# vars_needed <- c("tfs_months", "event_death", vars)
# df_sub_clinic <- df[complete.cases(df[, vars_needed]), vars_needed]
# 
# cat("Rows used:", nrow(df_sub_clinic), "\n")
# 
# # --- Build design matrix (one-hot encode factors) ---
# X <- model.matrix(~ . - 1, data = df_sub_clinic[, vars])  # removes intercept, encodes factors
# 
# # --- Response vector (Surv object) ---
# y <- Surv(df_sub_clinic$tfs_months, df_sub_clinic$event_death)
# 
# # --- Cross-validated LASSO Cox regression ---
# cvfit <- cv.glmnet(
#   X, y,
#   family = "cox",
#   alpha = 1,        # 1 = LASSO
#   nfolds = 10,      # 10-fold cross-validation
#   maxit = 100000    # in case dataset is large
# )
# 
# plot(cvfit)
# 
# cat("λ_min  :", cvfit$lambda.min,  "\n")
# cat("λ_1se  :", cvfit$lambda.1se,  "\n")
# 
# #Extract selected coefficients
# coef_min <- coef(cvfit, s = "lambda.min")
# sel_vars <- rownames(coef_min)[as.numeric(coef_min) != 0]
# 
# cat("\n Variables selected by LASSO (λ_min):\n")
# print(sel_vars)
# 
# #Tabulate the results
# lasso_clinic_results <- data.frame(
#   Variable = rownames(coef_min),
#   Coefficient = as.numeric(coef_min)
# )
# lasso_clinic_results <- subset(lasso_clinic_results, Coefficient != 0)
# print(lasso_clinic_results)
# 
# #Result is from 22 Vars, only 12 were selected with Coef !=0



# 2nd Lasso for Ramayan Variables

ram_vars <- c(
  "X202", "X205", "X208", "X211", "X213", "X216", "X219", "X222", "X225", "X228", "X231", "X234", "X237", "X239", 
  "X242", "X245", "X248", "X251", "X254", "X257", "X260", "X262", "X265", "X268", "X271", "X274", "X277", "X280", 
  "X282", "X285", "X288", "X291", "X294", "X297", "X300", "X302", "X305", "X308", "X311", "X314", "X317", "X320", 
  "X322", "X325", "X328", "X331", "X334", "X337", "X339", "X342", "X345", "X348", "X351", "X354", "X357", "X359", 
  "X362", "X365", "X368", "X371", "X374", "X376", "X379", "X382", "X385", "X388", "X390", "X393", "X396", "X399",
  "X402", "X405", "X407", "X410", "X413", "X416", "X419", "X421", "X424", "X427", "X430", "X433", "X436", "X438", 
  "X441", "X444", "X447", "X450", "X452", "X455", "X458", "X461", "X464", "X466", "X469", "X472", "X475", "X478", 
  "X480", "X483", "X486", "X489", "X492", "X494", "X497", "X500", "X503", "X505", "X508", "X511", "X514", "X517", 
  "X519", "X522", "X525", "X528", "X530", "X533", "X536", "X539", "X542", "X544", "X547", "X550", "X553", "X555", 
  "X558", "X561", "X564", "X566", "X569", "X572", "X575", "X578", "X580", "X583", "X586", "X589", "X591", "X594", 
  "X597", "X600", "X602", "X605", "X608", "X611", "X613", "X616", "X619", "X622", "X624", "X627", "X630", "X633", 
  "X635", "X638", "X641", "X644", "X646", "X649", "X652", "X654", "X657", "X660", "X663", "X665", "X668", "X671", 
  "X674", "X676", "X679", "X682", "X684", "X687", "X690", "X693", "X695", "X698", "X701", "X704", "X706", "X709", 
  "X712", "X714", "X717", "X720", "X723", "X725", "X728", "X731", "X733", "X736", "X739", "X742", "X744", "X747", 
  "X750", "X752", "X755", "X758", "X760", "X763", "X766", "X769", "X771", "X774", "X777", "X779", "X782", "X785", 
  "X787", "X790", "X793", "X795", "X798", "X801", "X804", "X806", "X809", "X812", "X814", "X817", "X820", "X822", 
  "X825", "X828", "X830", "X833", "X836", "X838", "X841", "X844", "X846", "X849", "X852", "X854", "X857", "X860", 
  "X862", "X865", "X868", "X870", "X873", "X876", "X878", "X881", "X884", "X886", "X889", "X892", "X894", "X897", 
  "X900", "X902", "X905", "X908", "X910", "X913", "X916", "X918", "X921", "X924", "X926", "X929", "X931", "X934", 
  "X937", "X939", "X942", "X945", "X947", "X950", "X953", "X955", "X958", "X960", "X963", "X966", "X968", "X971", 
  "X974", "X976", "X979", "X982", "X984", "X987", "X989", "X992", "X995", "X997", "X1000", "X1003", "X1005", "X1008",
  "X1010", "X1013", "X1016", "X1018", "X1021", "X1023", "X1026", "X1029", "X1031", "X1034", "X1037", "X1039", "X1042",
  "X1044", "X1047", "X1050", "X1052", "X1055", "X1057", "X1060", "X1063", "X1065", "X1068", "X1070", "X1073", "X1076",
  "X1078", "X1081", "X1083", "X1086", "X1089", "X1091", "X1094", "X1096", "X1099", "X1102", "X1104", "X1107", "X1109", 
  "X1112", "X1114", "X1117", "X1120", "X1122", "X1125", "X1127", "X1130", "X1133", "X1135", "X1138", "X1140", "X1143", 
  "X1145", "X1148", "X1151", "X1153", "X1156", "X1158", "X1161", "X1163", "X1166", "X1169", "X1171", "X1174", "X1176", 
  "X1179", "X1181", "X1184", "X1187", "X1189", "X1192", "X1194", "X1197", "X1199", "X1202", "X1205", "X1207", "X1210",
  "X1212", "X1215", "X1217", "X1220", "X1222", "X1225", "X1227", "X1230", "X1233", "X1235", "X1238", "X1240", "X1243", 
  "X1245", "X1248", "X1250", "X1253", "X1255", "X1258", "X1261", "X1263", "X1266", "X1268", "X1271", "X1273", "X1276", 
  "X1278", "X1281", "X1283", "X1286", "X1288", "X1291", "X1294", "X1296", "X1299", "X1301", "X1304", "X1306", "X1309", 
  "X1311", "X1314", "X1316", "X1319", "X1321", "X1324", "X1326", "X1329", "X1331", "X1334", "X1336", "X1339", "X1341",
  "X1344", "X1346", "X1349", "X1351", "X1354", "X1357", "X1359", "X1362", "X1364", "X1367", "X1369", "X1372", "X1374", 
  "X1377", "X1379", "X1382", "X1384", "X1387", "X1389", "X1392", "X1394", "X1397", "X1399", "X1402", "X1404", "X1407", 
  "X1409", "X1412", "X1414", "X1417", "X1419", "X1422", "X1424", "X1426", "X1429", "X1431", "X1434", "X1436", "X1439", 
  "X1441", "X1444", "X1446", "X1449", "X1451", "X1454", "X1456", "X1459", "X1461", "X1464", "X1466", "X1469", "X1471", 
  "X1474", "X1476", "X1479", "X1481", "X1484", "X1486", "X1488", "X1491", "X1493", "X1496", "X1498", "X1501", "X1503",
  "X1506", "X1508", "X1511", "X1513", "X1516", "X1518", "X1521", "X1523", "X1525", "X1528", "X1530", "X1533", "X1535",
  "X1538", "X1540", "X1543", "X1545", "X1548", "X1550", "X1552", "X1555", "X1557", "X1560", "X1562", "X1565", "X1567", 
  "X1570", "X1572", "X1574", "X1577", "X1579", "X1582", "X1584", "X1587", "X1589", "X1592", "X1594", "X1596", "X1599",
  "X1601", "X1604", "X1606", "X1609", "X1611", "X1613", "X1616", "X1618", "X1621", "X1623", "X1626", "X1628", "X1630", 
  "X1633", "X1635", "X1638", "X1640", "X1643", "X1645", "X1647", "X1650", "X1652", "X1655", "X1657", "X1659", "X1662",
  "X1664", "X1667", "X1669", "X1672", "X1674", "X1676", "X1679", "X1681", "X1684", "X1686", "X1688", "X1691", "X1693", 
  "X1696", "X1698", "X1700", "X1703", "X1705", "X1708", "X1710", "X1712", "X1715", "X1717", "X1720", "X1722", "X1724",
  "X1727", "X1729", "X1732", "X1734", "X1736", "X1739", "X1741", "X1744", "X1746", "X1748", "X1751", "X1753", "X1756", 
  "X1758", "X1760", "X1763", "X1765", "X1768", "X1770", "X1772", "X1775", "X1777", "X1779", "X1782", "X1784", "X1787", 
  "X1789", "X1791", "X1794", "X1796", "X1798", "X1801", "X1803", "X1806", "X1808", "X1810", "X1813", "X1815", "X1817", 
  "X1820", "X1822", "X1825", "X1827", "X1829", "X1832", "X1834", "X1836", "X1839", "X1841", "X1843", "X1846", "X1848", 
  "X1851", "X1853", "X1855", "X1858", "X1860", "X1862", "X1865", "X1867", "X1869", "X1872", "X1874", "X1876", "X1879", 
  "X1881", "X1883", "X1886", "X1888", "X1891", "X1893", "X1895", "X1898", "X1900", "X1902", "X1905", "X1907", "X1909",
  "X1912", "X1914", "X1916", "X1919", "X1921", "X1923", "X1926", "X1928", "X1930", "X1933", "X1935", "X1937", "X1940", 
  "X1942", "X1944", "X1947", "X1949", "X1951", "X1954", "X1956", "X1958", "X1961", "X1963", "X1965", "X1968", "X1970",
  "X1972", "X1975", "X1977", "X1979", "X1982", "X1984", "X1986", "X1989", "X1991", "X1993", "X1995", "X1998", "X2000",
  "X2002", "X2005", "X2007", "X2009", "X2012", "X2014", "X2016", "X2019", "X2021", "X2023", "X2026", "X2028", "X2030",
  "X2032", "X2035", "X2037", "X2039", "X2042", "X2044", "X2046", "X2049", "X2051", "X2053", "X2055", "X2058", "X2060", 
  "X2062", "X2065", "X2067", "X2069", "X2072", "X2074", "X2076", "X2078", "X2081", "X2083", "X2085", "X2088", "X2090",
  "X2092", "X2094", "X2097", "X2099", "X2101", "X2104", "X2106", "X2108", "X2110", "X2113", "X2115", "X2117", "X2120",
  "X2122", "X2124", "X2126", "X2129", "X2131", "X2133", "X2135", "X2138", "X2140", "X2142", "X2145", "X2147", "X2149",
  "X2151", "X2154", "X2156", "X2158", "X2160", "X2163", "X2165", "X2167", "X2170", "X2172", "X2174", "X2176", "X2179",
  "X2181", "X2183", "X2185", "X2188", "X2190", "X2192", "X2194", "X2197", "X2199", "X2201", "X2203", "X2206", "X2208",
  "X2210", "X2212", "X2215", "X2217", "X2219", "X2221", "X2224", "X2226", "X2228", "X2230", "X2233", "X2235", "X2237",
  "X2239", "X2242", "X2244", "X2246", "X2248", "X2251", "X2253", "X2255", "X2257", "X2260", "X2262", "X2264", "X2266",
  "X2269", "X2271", "X2273", "X2275", "X2278", "X2280", "X2282", "X2284", "X2286", "X2289", "X2291", "X2293", "X2295",
  "X2298", "X2300", "X2302", "X2304", "X2306", "X2309", "X2311", "X2313", "X2315", "X2318", "X2320", "X2322", "X2324", 
  "X2326", "X2329", "X2331", "X2333", "X2335", "X2338", "X2340", "X2342", "X2344", "X2346", "X2349", "X2351", "X2353", 
  "X2355", "X2358", "X2360", "X2362", "X2364", "X2366", "X2369", "X2371", "X2373", "X2375", "X2377", "X2380", "X2382", 
  "X2384", "X2386", "X2388", "X2391", "X2393", "X2395", "X2397", "X2399", "X2402", "X2404", "X2406", "X2408", "X2410", 
  "X2413", "X2415", "X2417", "X2419", "X2421", "X2424", "X2426", "X2428", "X2430", "X2432", "X2434", "X2437", "X2439",
  "X2441", "X2443", "X2445", "X2448", "X2450", "X2452", "X2454", "X2456", "X2459", "X2461", "X2463", "X2465", "X2467", 
  "X2469", "X2472", "X2474", "X2476", "X2478", "X2480", "X2482", "X2485", "X2487", "X2489", "X2491", "X2493", "X2496",
  "X2498", "X2500", "X2502", "X2504", "X2506", "X2509", "X2511", "X2513", "X2515", "X2517", "X2519", "X2522", "X2524", 
  "X2526", "X2528", "X2530", "X2532", "X2534", "X2537", "X2539", "X2541", "X2543", "X2545", "X2547", "X2550", "X2552",
  "X2554", "X2556", "X2558", "X2560", "X2563", "X2565", "X2567", "X2569", "X2571", "X2573", "X2575", "X2578", "X2580", 
  "X2582", "X2584", "X2586", "X2588", "X2590", "X2593", "X2595", "X2597", "X2599", "X2601", "X2603", "X2605", "X2608", 
  "X2610", "X2612", "X2614", "X2616", "X2618", "X2620", "X2623", "X2625", "X2627", "X2629", "X2631", "X2633", "X2635", 
  "X2637", "X2640", "X2642", "X2644", "X2646", "X2648", "X2650", "X2652", "X2654", "X2657", "X2659", "X2661", "X2663", 
  "X2665", "X2667", "X2669", "X2671", "X2674", "X2676", "X2678", "X2680", "X2682", "X2684", "X2686", "X2688", "X2691", 
  "X2693", "X2695", "X2697", "X2699", "X2701", "X2703", "X2705", "X2707", "X2710", "X2712", "X2714", "X2716", "X2718", 
  "X2720", "X2722", "X2724", "X2726", "X2729", "X2731", "X2733", "X2735", "X2737", "X2739", "X2741", "X2743", "X2745", 
  "X2747", "X2750", "X2752", "X2754", "X2756", "X2758", "X2760", "X2762", "X2764", "X2766", "X2768", "X2771", "X2773", 
  "X2775", "X2777", "X2779", "X2781", "X2783", "X2785", "X2787", "X2789", "X2791", "X2793", "X2796", "X2798", "X2800", 
  "X2802", "X2804", "X2806", "X2808", "X2810", "X2812", "X2814", "X2816", "X2818", "X2821", "X2823", "X2825", "X2827", 
  "X2829", "X2831", "X2833", "X2835", "X2837", "X2839", "X2841", "X2843", "X2845", "X2848", "X2850", "X2852", "X2854", 
  "X2856", "X2858", "X2860", "X2862", "X2864", "X2866", "X2868", "X2870", "X2872", "X2874", "X2877", "X2879", "X2881", 
  "X2883", "X2885", "X2887", "X2889", "X2891", "X2893", "X2895", "X2897", "X2899", "X2901", "X2903", "X2905", "X2907", 
  "X2909", "X2912", "X2914", "X2916", "X2918", "X2920", "X2922", "X2924", "X2926", "X2928", "X2930", "X2932", "X2934", 
  "X2936", "X2938", "X2940", "X2942", "X2944", "X2946", "X2948", "X2950", "X2952", "X2955", "X2957", "X2959", "X2961", 
  "X2963", "X2965", "X2967", "X2969", "X2971", "X2973", "X2975", "X2977", "X2979", "X2981", "X2983", "X2985", "X2987", 
  "X2989", "X2991", "X2993", "X2995", "X2997", "X2999", "X3001", "X3003", "X3005", "X3007", "X3009", "X3012", "X3014", 
  "X3016", "X3018", "X3020", "X3022", "X3024", "X3026", "X3028", "X3030", "X3032", "X3034", "X3036", "X3038", "X3040", 
  "X3042", "X3044", "X3046", "X3048", "X3050", "X3052", "X3054", "X3056", "X3058", "X3060", "X3062", "X3064", "X3066", 
  "X3068", "X3070", "X3072", "X3074", "X3076", "X3078", "X3080", "X3082", "X3084", "X3086", "X3088", "X3090", "X3092", 
  "X3094", "X3096", "X3098", "X3100", "X3102", "X3104", "X3106", "X3108", "X3110", "X3112", "X3114", "X3116", "X3118", 
  "X3120", "X3122", "X3124", "X3126", "X3128", "X3130", "X3132", "X3134", "X3136", "X3138", "X3140", "X3142", "X3144", 
  "X3146", "X3148"
)




# df_ram_filtered <- df[complete.cases(df[, ram_vars]), ]
# 
# # Print how many rows remain and how many were removed
# cat("Rows remaining :", nrow(df_ram_filtered),"\n")
# 
# set.seed(123)  # for reproducibility
# 
# # Define event and predictors 
# y <- Surv(df$tfs_months, df$event_death)
# vars <- ram_vars
# 
# # Keep only rows with no NA in these variables and survival columns
# vars_needed <- c("tfs_months", "event_death", vars)
# df_sub_rama <- df[complete.cases(df[, vars_needed]), vars_needed]
# length(vars)
# dim(df_sub_rama)
# cat("Rows used:", nrow(df_sub_rama), "\n")
# 
# 
# # --- Build design matrix (one-hot encode factors) ---
# X_rama <- model.matrix(~ . - 1, data = df_sub_rama[, vars])  # removes intercept, encodes factors
# 
# # --- Response vector (Surv object) ---
# y <- Surv(df_sub_rama$tfs_months, df_sub_rama$event_death)
# 
# # --- Cross-validated LASSO Cox regression ---
# cvfit_rama <- cv.glmnet(
#   X_rama, y,
#   family = "cox",
#   alpha = 1,        # 1 = LASSO
#   nfolds = 10,      # 10-fold cross-validation
#   maxit = 100000    # in case dataset is large
# )
# 
# plot(cvfit_rama)
# 
# cat("λ_min  :", cvfit_rama$lambda.min,  "\n")
# cat("λ_1se  :", cvfit_rama$lambda.1se,  "\n")
# 
# #Extract selected coefficients
# coef_min <- coef(cvfit_rama, s = "lambda.min")
# sel_vars <- rownames(coef_min)[as.numeric(coef_min) != 0]
# 
# cat("\n Variables selected by LASSO (λ_min):\n")
# print(sel_vars)
# 
# #Tabulate the results
# lasso_rama_results <- data.frame(
#   Variable = rownames(coef_min),
#   Coefficient = as.numeric(coef_min)
# )
# lasso_rama_results <- subset(lasso_rama_results, Coefficient != 0)
# print(lasso_rama_results)
# 
# #Result is that from 1231 Ramayan Vars, only 19 have coef !=0


# 3st do Lasso for Lab variables
# From the Lab Data ,hba1c, iga, igg, igm, igg_1, igg_2, igg_3, igg_4 wont be used in the regression.
lab_vars <- c(
  "hb", "plt", "wbc", "na", "k", "cl", "krea", "bun", "bili", "protein", "alb", "che", "ap",
  "asat", "alat", "ggt", "ldh", "aamy", "pamy", "lip", "vwf", "fgen", "inr", "ba", "trig", 
  "chol", "gluc", "hdl", "ldl","nh3", "crp", "il6", "pct", "lbp", "copeptinadh", 
  "bnp", "elf", "timp1", "p3np", "ha"
)

lab_all_vars<- c(
  "hb", "plt", "wbc", "na", "k", "cl", "krea", "bun", "bili", "protein", "alb", "che", "ap",
  "asat", "alat", "ggt", "ldh", "aamy", "pamy", "lip", "vwf", "fgen", "inr", "ba", "trig", 
  "chol", "gluc", "hdl", "ldl", "hba1c", "nh3", "crp", "il6", "pct", "lbp", "copeptinadh", 
  "bnp", "elf", "timp1", "p3np", "ha", "iga", "igg", "igm", "igg_1", "igg_2", "igg_3", "igg_4"
)

# df_lab<- df[lab_vars]
# 
# df_lab_filtered <- df[complete.cases(df[, lab_vars]), ]
# 
# # Print how many rows remain and how many were removed
# cat("Rows remaining :", nrow(df_lab_filtered),"\n")
# 
# 
# set.seed(123)  # for reproducibility
# 
# # Define event and predictors 
# y <- Surv(df$tfs_months, df$event_death)
# vars <- lab_vars
# 
# # Keep only rows with no NA in these variables and survival columns
# vars_needed <- c("tfs_months", "event_death", vars)
# df_sub_lab <- df[complete.cases(df[, vars_needed]), vars_needed]
# 
# cat("Rows used:", nrow(df_sub_lab), "\n")
# 
# # --- Build design matrix (one-hot encode factors) ---
# X <- model.matrix(~ . - 1, data = df_sub_lab[, vars])  # removes intercept, encodes factors
# 
# # --- Response vector (Surv object) ---
# y <- Surv(df_sub_lab$tfs_months, df_sub_lab$event_death)
# 
# # --- Cross-validated LASSO Cox regression ---
# cvfit <- cv.glmnet(
#   X, y,
#   family = "cox",
#   alpha = 1,        # 1 = LASSO
#   nfolds = 10,      # 10-fold cross-validation
#   maxit = 100000    # in case dataset is large
# )
# 
# plot(cvfit)
# 
# cat("λ_min  :", cvfit$lambda.min,  "\n")
# cat("λ_1se  :", cvfit$lambda.1se,  "\n")
# 
# #Extract selected coefficients
# coef_min <- coef(cvfit, s = "lambda.min")
# sel_vars <- rownames(coef_min)[as.numeric(coef_min) != 0]
# 
# cat("\n Variables selected by LASSO (λ_min):\n")
# print(sel_vars)
# 
# #Tabulate the results
# lasso_lab_results <- data.frame(
#   Variable = rownames(coef_min),
#   Coefficient = as.numeric(coef_min)
# )
# lasso_lab_results <- subset(lasso_lab_results, Coefficient != 0)
# print(lasso_lab_results)
# 
# #Result is from 47 Vars, only 4 were selected with Coef !=0

# 4th do Lasso for Bile Acid variables
# From the BA Data biobank_nr, biobank_kerstin, tmca, tcasulfate, tudca, omca,casulfate
#amca, bmca, udca, tlca,lca,  wont be used in the regression.
ba_vars <- c(
    "tca", "gca", "gudca", "tcdca", "tdca", "ca", 
  "gcdca", "gdca", "cdca", "glca","dca", "c4", "totalba", 
  "totalba_muM", "kim1", "ngal", "cystatinc"
)


#df_ba_filtered <- df[complete.cases(df[, ba_vars]), ]
# 
# # Print how many rows remain and how many were removed
# cat("Rows remaining :", nrow(df_ba_filtered),"\n")
# 
# 
# set.seed(123)  # for reproducibility
# 
# # Define event and predictors 
# y <- Surv(df$tfs_months, df$event_death)
# vars <- ba_vars
# 
# # Keep only rows with no NA in these variables and survival columns
# vars_needed <- c("tfs_months", "event_death", vars)
# df_sub_ba <- df[complete.cases(df[, vars_needed]), vars_needed]
# 
# cat("Rows used:", nrow(df_sub_ba), "\n")
# 
# # Build design matrix (one-hot encode factors) 
# X <- model.matrix(~ . - 1, data = df_sub_ba[, vars])  # removes intercept, encodes factors
# 
# #  Response vector (Surv object)
# y <- Surv(df_sub_ba$tfs_months, df_sub_ba$event_death)
# 
# # Cross-validated LASSO Cox regression 
# cvfit <- cv.glmnet(
#   X, y,
#   family = "cox",
#   alpha = 1,        # 1 = LASSO
#   nfolds = 10,      # 10-fold cross-validation
#   maxit = 100000    # in case dataset is large
# )
# 
# plot(cvfit)
# 
# cat("λ_min  :", cvfit$lambda.min,  "\n")
# cat("λ_1se  :", cvfit$lambda.1se,  "\n")
# 
# #Extract selected coefficients
# coef_min <- coef(cvfit, s = "lambda.min")
# sel_vars <- rownames(coef_min)[as.numeric(coef_min) != 0]
# 
# cat("\n Variables selected by LASSO (λ_min):\n")
# print(sel_vars)
# 
# #Tabulate the results
# lasso_ba_results <- data.frame(
#   Variable = rownames(coef_min),
#   Coefficient = as.numeric(coef_min)
# )
# lasso_ba_results <- subset(lasso_ba_results, Coefficient != 0)
# print(lasso_ba_results)
# 
# #Result is from 21 Vars, only 3 were selected with Coef !=0




#-------------------------------------------------------------------------------
# Descriptive survival analysis (Kaplan–Meier, summaries, group comparison)

library(survival)  # Loads core survival-analysis functions: Surv(), survfit(), survdiff(), etc.


## -------------------------
## 1) Create survival object
## -------------------------

# Creates a survival response object combining:
# - time: follow-up time (months)
# - event: 1=event occurred (death here), 0=censored
# This is the standard input format for Kaplan–Meier and Cox models.
surv_obj_desc <- Surv(
  time  = df$tfs_months,
  event = df$event_death
)


## -------------------------
## 2) Overall Kaplan–Meier
## -------------------------

# Fits an overall Kaplan–Meier curve (no groups) because the formula is "~ 1".
# This estimates S(t) = P(T > t) nonparametrically, accounting for censoring.
km_fit_desc <- survfit(
  surv_obj_desc ~ 1,
  data = df
)

# Plots the KM survival curve:
# - conf.int=TRUE adds confidence bands
# - mark.time=TRUE marks censoring times on the curve
plot(
  km_fit_desc,
  xlab = "Time (months)",
  ylab = "Survival probability",
  main = "Kaplan–Meier Survival Curve",
  conf.int = TRUE,
  mark.time = TRUE
)

# Extracts detailed KM output (time points, survival estimates, CI, # at risk, # events).
km_summary_desc <- summary(km_fit_desc)

# Shows the "table" component with key overall summaries:
# typically includes n, events, median survival (if reached), and CI for median.
km_summary_desc$table


## -------------------------
## 3) Follow-up time summaries
## -------------------------

# Quick numeric summary of follow-up time (regardless of event/censoring):
# min, 1st quartile, median, mean, 3rd quartile, max
surv_time_summary_desc <- summary(df$tfs_months)
surv_time_summary_desc

# Summarizes follow-up time separately by event status:
# event_death==1: observed event times
# event_death==0: censored times (follow-up until last contact)
# Useful to see if censoring happens earlier/later than events.
surv_time_by_event_desc <- with(
  df,
  tapply(tfs_months, event_death, summary)
)
surv_time_by_event_desc


## -------------------------
## 4) Event counts / proportions
## -------------------------

# Counts how many are events vs censored.
event_table_desc <- table(df$event_death)
event_table_desc

# Converts counts into proportions (sums to 1).
# Helpful for reporting event rate / censoring rate.
event_prop_desc <- prop.table(event_table_desc)
event_prop_desc


## -------------------------
## 5) Kaplan–Meier stratified by sex + plot
## -------------------------

# KM curves by sex:
# Fits separate survival curves per sex group and estimates S(t) in each group.
km_fit_by_sex_desc <- survfit(
  Surv(tfs_months, event_death) ~ sex,
  data = df
)

# Plots both sex-specific curves on the same axes.
plot(
  km_fit_by_sex_desc,
  col = c("blue", "red"),
  lwd = 2,
  xlab = "Time (months)",
  ylab = "Survival probability",
  main = "Kaplan–Meier by Sex",
  conf.int = TRUE,     # adds CI bands for each group
  mark.time = TRUE     # marks censoring times
)

# Adds legend mapping colors to the sex factor levels.
legend(
  "bottomleft",
  legend = levels(df$sex),
  col = c("blue", "red"),
  lwd = 2
)


## -------------------------
## 6) Log-rank test (group comparison)
## -------------------------

# Performs a log-rank test comparing survival curves between sexes.
# Null hypothesis: the survival functions are identical across groups.
# Output includes observed vs expected events and a chi-square test statistic with p-value.
logrank_sex_desc <- survdiff(
  Surv(tfs_months, event_death) ~ sex,
  data = df
)
logrank_sex_desc


## -------------------------
## 7) Histograms of times by event status
## -------------------------

# Histogram of observed event times (death times).
# probability=TRUE plots density (area under bars sums to 1).
# Semi-transparent red so we can overlay with the censored distribution.
hist(
  df$tfs_months[df$event_death == 1],
  breaks = 30,
  col = rgb(1, 0, 0, 0.5),
  probability = TRUE,
  main = "Survival Times by Event Status",
  xlab = "Time (months)"
)

# Histogram of censoring times (follow-up until censoring), overlaid on the same plot.
# add=TRUE overlays the second histogram onto the first.
hist(
  df$tfs_months[df$event_death == 0],
  breaks = 30,
  probability = TRUE,
  col = rgb(0, 0, 1, 0.5),
  add = TRUE
)

# Legend explaining the color coding for the two distributions.
legend(
  "topright",
  legend = c("Event", "Censored"),
  fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))
)


#-------------------------------------------------------------------------------
# Double CV to calc BS and IPA

#simple imputation for collumns with Nas, if continous then use median and if 
#categorical then mode





time_cols  <- c("hea_months","reche_months", "aki_months","hrsaki_months","inf_months"
                ,"firdec_cacld_months","furdec_cacld_months","furdec_dacld_months"
                ,"aclf_months","lrd_months","aclf_lrd_months","cumulative_months"
                ,"reg_dacld_strict_months" ,"tfs_months")
event_cols <- c("hea_01","reche_01", "aki_01", "hrsaki_01","inf_01","firdec_cacld_01",
                "furdec_cacld_01" ,"furdec_dacld_01","aclf_01", "lrd_013","aclf_lrd_01","cumulative_01"
                ,"reg_dacld_strict_01", "tfs_event_012")
surv_cols  <- c(time_cols, event_cols)


# Identify columns with >14% missing (excluding survival cols)

na_fraction <- sapply(df, function(x) mean(is.na(x)))

cols_remove <- names(na_fraction[
  na_fraction > 0.14 & !(names(na_fraction) %in% surv_cols)
])

cols_keep <- setdiff(names(df), cols_remove)

cat("Columns removed (>14% NA):\n")
print(cols_remove)


# Keep only surviving columns

df_clean <- df[, cols_keep, drop = FALSE]

#Impute all remaining variables EXCEPT the survival columns

impute_simple <- function(df, surv_cols) {
  
  for (col in names(df)) {
    
    # skip survival columns entirely
    if (col %in% surv_cols) next
    
    x <- df[[col]]
    
    # numeric → median
    if (is.numeric(x)) {
      med <- median(x, na.rm = TRUE)
      x[is.na(x)] <- med
      df[[col]] <- x
    }
    
    # factor/character → mode
    else {
      tab <- table(x)
      mode_val <- names(tab)[which.max(tab)]
      x[is.na(x)] <- mode_val
      df[[col]] <- x
    }
  }
  
  return(df)
}



# Apply imputation
df_imputed <- impute_simple(df_clean, surv_cols)

cat("\nImputation complete.\n")
cat("Final dataset dimensions: ")
print(dim(df_imputed))



#create combos of vars type , eg ba+ram or lab+ba+clinical
make_var_group_combos <- function(var_groups, sizes = 1:length(var_groups)) {
  group_names <- names(var_groups)
  combos <- list()
  id <- 1
  for (k in sizes) {
    cmb <- utils::combn(group_names, k, simplify = FALSE)
    for (g in cmb) {
      combos[[id]] <- g
      id <- id + 1
    }
  }
  combos
}


# Main function that conducts the inner(to find the optimal lambda) and outer(
# to calc the BS and IPA) cv.
run_lasso_over_combos_doublecv <- function(
    df,
    var_groups,
    time_col    = "tfs_months",
    event_col   = "event_death",
    K_outer     = 10,
    K_inner     = 10,
    combo_sizes = 1:length(var_groups),
    maxit       = 100000,
    min_rows    = 30,
    min_events  = 10,
    times       = seq(12, 60, 12)
) {

  combos <- make_var_group_combos(var_groups, sizes = combo_sizes)
  summary_rows <- list()
  details_rows <- list()
  combo_id <- 1L

  df[[time_col]] <- ifelse(df[[time_col]] <= 0, 1e-6, df[[time_col]])

  for (groups_in in combos) {

    vars <- unique(unlist(var_groups[groups_in], use.names = FALSE))
    combo_label <- paste(groups_in, collapse = "+")

    vars_needed <- c(time_col, event_col, vars)
    missing_cols <- setdiff(vars_needed, colnames(df))

    if (length(missing_cols) > 0) {
      summary_rows[[length(summary_rows)+1]] <- data.frame(
        combo_id = combo_id,
        groups = combo_label,
        n_vars = length(vars),
        rows_used = NA,
        events = NA,
        IPA_mean = NA,
        IPA_sd = NA,
        BS_mean = NA,
        BS_sd = NA,
        status = paste("missing:", paste(missing_cols, collapse=",")),
        stringsAsFactors = FALSE
      )
      combo_id <- combo_id + 1L
      next
    }

    dsub <- df[complete.cases(df[, vars_needed]), vars_needed, drop = FALSE]
    dsub <- dsub[order(rownames(dsub)), , drop = FALSE]
    n_rows <- nrow(dsub)
    n_events <- sum(dsub[[event_col]] == 1)

    if (n_rows < min_rows || n_events < min_events) {
      summary_rows[[length(summary_rows)+1]] <- data.frame(
        combo_id = combo_id,
        groups = combo_label,
        n_vars = length(vars),
        rows_used = n_rows,
        events = n_events,
        IPA_mean = NA,
        IPA_sd = NA,
        BS_mean = NA,
        BS_sd = NA,
        status = "insufficient rows/events",
        stringsAsFactors = FALSE
      )
      combo_id <- combo_id + 1L
      next
    }

    X_all <- model.matrix(~ . - 1, data = dsub[, vars, drop = FALSE])
    y_all <- Surv(dsub[[time_col]], dsub[[event_col]])

    set.seed(123)
    e <- as.integer(dsub[[event_col]])
    fold_outer <- integer(n_rows)
    idx1 <- which(e == 1)
    idx0 <- which(e == 0)
    fold_outer[idx1] <- sample(rep(1:K_outer, length.out = length(idx1)))
    fold_outer[idx0] <- sample(rep(1:K_outer, length.out = length(idx0)))

    ipa_vec <- rep(NA_real_, K_outer)
    bs_vec <- rep(NA_real_, K_outer)
    bs_null_vec <- rep(NA_real_, K_outer)
    selected_vars_per_fold <- vector("list", K_outer)
    lambda_min_vec <- rep(NA_real_, K_outer)
    

    surv_formula_score <- Surv(time, event) ~ 1

    for (o in 1:K_outer) {

      train_idx <- which(fold_outer != o)
      test_idx <- which(fold_outer == o)

      df_train <- dsub[train_idx, , drop = FALSE]
      df_test <- dsub[test_idx, , drop = FALSE]
      X_train <- X_all[train_idx, , drop = FALSE]
      X_test <- X_all[test_idx, , drop = FALSE]

      y_train <- Surv(df_train[[time_col]], df_train[[event_col]])

      cvfit <- cv.glmnet(
        x = X_train,
        y = y_train,
        family = "cox",
        alpha = 1,
        nfolds = K_inner,
        maxit = maxit
      )
      
      lambda_min_vec[o] <- cvfit$lambda.min

      coef_min <- coef(cvfit, s = "lambda.min")
      sel_idx <- as.numeric(coef_min) != 0
      sel_names <- rownames(coef_min)[sel_idx]
      sel_names <- setdiff(sel_names, "(Intercept)")

      selected_vars_per_fold[[o]] <- sel_names

      if (length(sel_names) == 0) {
        ipa_vec[o] <- NA
        bs_vec[o] <- NA
        bs_null_vec[o] <- NA
        next
      }

      X_train_sel <- X_train[, sel_names, drop = FALSE]
      X_test_sel <- X_test[, sel_names, drop = FALSE]

      df_train_cox <- data.frame(
        time = df_train[[time_col]],
        event = df_train[[event_col]],
        X_train_sel,
        check.names = FALSE
      )

      df_test_cox <- data.frame(
        time = df_test[[time_col]],
        event = df_test[[event_col]],
        X_test_sel,
        check.names = FALSE
      )

      safe_names <- paste0("`", colnames(X_train_sel), "`")

      form_cox <- as.formula(
        paste0("Surv(time, event) ~ ", paste(safe_names, collapse = " + "))
      )

      fit_outer <- try(
        coxph(form_cox, data = df_train_cox, x = TRUE),
        silent = TRUE
      )

      if (inherits(fit_outer, "try-error")) {
        ipa_vec[o] <- NA
        bs_vec[o] <- NA
        bs_null_vec[o] <- NA
        next
      }

      if (any(is.na(fit_outer$coefficients))) {

        bad_vars <- names(fit_outer$coefficients)[is.na(fit_outer$coefficients)]
        sel_names_clean <- setdiff(sel_names, bad_vars)

        if (length(sel_names_clean) == 0) {
          ipa_vec[o] <- NA
          bs_vec[o] <- NA
          bs_null_vec[o] <- NA
          next
        }

        X_train_sel <- X_train[, sel_names_clean, drop = FALSE]
        X_test_sel <- X_test[, sel_names_clean, drop = FALSE]

        df_train_cox <- data.frame(
          time = df_train[[time_col]],
          event = df_train[[event_col]],
          X_train_sel,
          check.names = FALSE
        )

        df_test_cox <- data.frame(
          time = df_test[[time_col]],
          event = df_test[[event_col]],
          X_test_sel,
          check.names = FALSE
        )

        safe_names <- paste0("`", colnames(X_train_sel), "`")

        form_cox <- as.formula(
          paste0("Surv(time, event) ~ ", paste(safe_names, collapse = " + "))
        )

        fit_outer <- try(
          coxph(form_cox, data = df_train_cox, x = TRUE),
          silent = TRUE
        )

        if (inherits(fit_outer, "try-error")) {
          ipa_vec[o] <- NA
          bs_vec[o] <- NA
          bs_null_vec[o] <- NA
          next
        }
      }

      max_test_time <- max(df_test_cox$time, na.rm = TRUE)
      times_fold <- times[times <= max_test_time]
      
      if (length(times_fold) == 0) {
        ipa_vec[o] <- NA
        bs_vec[o] <- NA
        bs_null_vec[o] <- NA
        next
      }
      

      score_obj <- try(
        Score(
          object = list("lasso" = fit_outer),
          formula = surv_formula_score,
          data = df_test_cox,
          metrics = "Brier",
          times = times_fold,
          null.model = TRUE,
          summary = "IPA"
        ),
        silent = TRUE
      )

      if (inherits(score_obj, "try-error")) {
        ipa_vec[o] <- NA
        bs_vec[o] <- NA
        bs_null_vec[o] <- NA
        next
      }

      BS_model_vec <- score_obj$Brier$score$Brier[
        score_obj$Brier$score$model == "lasso"
      ]
      
      BS_KM_vec <- score_obj$Brier$score$Brier[
        score_obj$Brier$score$model == "Null model"
      ]
      
      bs_vec[o]      <- mean(BS_model_vec, na.rm = TRUE)
      bs_null_vec[o] <- mean(BS_KM_vec,    na.rm = TRUE)
      
      ipa_vec[o] <- 1 - (bs_vec[o] / bs_null_vec[o])
      ### Debug output ###
      cat("\n*** COMBO", combo_label, "FOLD", o, "***")
      cat("\ntimes_fold:", paste(times_fold, collapse=", "))
      cat("\nBS_model:", bs_vec[o])
      cat("\nBS_KM:", bs_null_vec[o])
      cat("\nIPA_fold:", ipa_vec[o], "\n")

    }  # END outer fold loop

    summary_rows[[length(summary_rows)+1]] <- data.frame(
      combo_id = combo_id,
      groups = combo_label,
      n_vars = length(vars),
      rows_used = n_rows,
      events = n_events,
      IPA_mean = mean(ipa_vec, na.rm = TRUE),
      IPA_sd = sd(ipa_vec, na.rm = TRUE),
      BS_mean = mean(bs_vec, na.rm = TRUE),
      BS_sd = sd(bs_vec, na.rm = TRUE),
      BS_mean_null = mean(bs_null_vec, na.rm = TRUE), 
      BS_sd_null = sd(bs_null_vec, na.rm = TRUE),
      status = "ok",
      stringsAsFactors = FALSE
    )

    details_rows[[as.character(combo_id)]] <- list(
      groups = groups_in,
      vars = vars,
      selected_vars_per_fold = selected_vars_per_fold,
      lambda_min = lambda_min_vec,
      IPA = ipa_vec,
      BS = bs_vec,
      BS_null= bs_null_vec
    )

    combo_id <- combo_id + 1L
  }  # END combo loop

  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL

  list(summary = summary_df, details = details_rows)
}


var_groups <- list(
  lab      = lab_vars,
  ba       = ba_vars,
  ram      = ram_vars,
  clinical = clinical_vars
)

df$event_death <- ifelse(df$tfs_event_012 == 2L, 1L, 0L)
df$tfs_months  <- ifelse(df$tfs_months <= 0, 1e-6, df$tfs_months)

#imputed but both uses the same vars. just that the imputed has imputation on NAs
results_doublecv_imputed <- run_lasso_over_combos_doublecv(
  df          = df_imputed,
  var_groups  = var_groups,
  combo_sizes = 1:4,
  K_outer     = 10,
  K_inner     = 10,
  min_rows    = 150,
  min_events  = 10,
  times       = seq(12, 60, 12)
)
results_doublecv_imputed$summary

#nonimputed df

# results_doublecv <- run_lasso_over_combos_doublecv(
#   df          = df,
#   var_groups  = var_groups,
#   combo_sizes = 1:4,
#   K_outer     = 10,
#   K_inner     = 10,
#   min_rows    = 150,
#   min_events  = 10,
#   times       = seq(12, 60, 12)
# )
# results_doublecv$summary





#Heatmap of NAs with raman vars taken out since theres no Nas
# all_collumns<- colnames(df)
# all_cols_without_ram<- setdiff(all_collumns,ram_vars)
# df_sub_col<- df[,all_cols_without_ram, drop= FALSE]
# na_matrix<- is.na(df_sub_col)
# na_numeric<- na_matrix*1
# heatmap(na_numeric)
# 
# all_collumns<- colnames(df_imputed)
# all_cols_without_ram<- setdiff(all_collumns,ram_vars)
# df_sub_col<- df[,all_cols_without_ram, drop= FALSE]
# na_matrix<- is.na(df_sub_col)
# na_numeric<- na_matrix*1
# heatmap(na_numeric)




#-----------------------------------------------------
#SRF
#-----------------------------------------------------


doublecv_srf_brier <- function(
    df,
    vars,
    time_col  = "tfs_months",
    event_col = "event_death",
    K_outer   = 10,
    K_inner   = 5,
    times     = seq(12, 60, 12),
    
    # inner tuning grid
    grid = expand.grid(
      mtry     = c(2, 4, 8),
      nodesize = c(5, 10, 15),
      nsplit   = c(10),
      ntree    = c(500),
      stringsAsFactors = FALSE
    ),
    
    seed      = 123,
    verbose   = TRUE
) {
  
  # -----------------------------
  # Complete-case subset
  # -----------------------------
  needed <- c(time_col, event_col, vars)
  dsub <- df[complete.cases(df[, needed]), needed, drop = FALSE]
 
  if (nrow(dsub) == 0) stop("No rows left after complete.cases().")
  
  dsub[[time_col]] <- ifelse(dsub[[time_col]] <= 0, 1e-6, dsub[[time_col]])
  dsub <- dsub[order(rownames(dsub)), , drop = FALSE]
  # SRF formula
  srf_formula <- as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(vars, collapse = " + "))
  )
  
  # -----------------------------
  # Create OUTER folds (stratified by event)
  # -----------------------------
  set.seed(seed)
  n <- nrow(dsub)
  e <- as.integer(dsub[[event_col]])
  
  folds_outer <- integer(n)
  idx1 <- which(e == 1)
  idx0 <- which(e == 0)
  folds_outer[idx1] <- sample(rep(1:K_outer, length.out = length(idx1)))
  folds_outer[idx0] <- sample(rep(1:K_outer, length.out = length(idx0)))
  
  # -----------------------------
  # Storage for OUTER results + VIMP
  # -----------------------------
  bs_outer      <- rep(NA_real_, K_outer)
  bs_null_outer <- rep(NA_real_, K_outer)
  ipa_outer     <- rep(NA_real_, K_outer)
  
  best_params_outer <- vector("list", K_outer)
  inner_grid_tables <- vector("list", K_outer)
  
  # store fold-wise VIMP vectors
  vimp_outer <- vector("list", K_outer)
  
  # ==========================================================
  # OUTER CV LOOP
  # ==========================================================
  for (o in 1:K_outer) {
    
    if (verbose) message("\n===== OUTER FOLD ", o, "/", K_outer, " =====")
    
    train_idx <- which(folds_outer != o)
    test_idx  <- which(folds_outer == o)
    
    dtrain_outer <- dsub[train_idx, , drop = FALSE]
    dtest_outer  <- dsub[test_idx,  , drop = FALSE]
    
    # basic event sanity
    if (sum(dtrain_outer[[event_col]] == 1) < 2 || sum(dtest_outer[[event_col]] == 1) < 1) {
      if (verbose) message("Outer fold ", o, ": too few events; skipping.")
      next
    }
    
    # restrict evaluation times to observed range in OUTER test fold
    max_test_time <- max(dtest_outer[[time_col]], na.rm = TRUE)
    times_outer <- times[times <= max_test_time]
    if (length(times_outer) == 0) {
      if (verbose) message("Outer fold ", o, ": no usable evaluation times; skipping.")
      next
    }
    
    # ----------------------------------------------------------
    # INNER folds (stratified) within OUTER training set
    # ----------------------------------------------------------
    set.seed(seed + o)
    ntr <- nrow(dtrain_outer)
    e_tr <- as.integer(dtrain_outer[[event_col]])
    
    folds_inner <- integer(ntr)
    id1 <- which(e_tr == 1)
    id0 <- which(e_tr == 0)
    folds_inner[id1] <- sample(rep(1:K_inner, length.out = length(id1)))
    folds_inner[id0] <- sample(rep(1:K_inner, length.out = length(id0)))
    
    # ----------------------------------------------------------
    # INNER CV: tune hyperparameters by mean Brier score
    # ----------------------------------------------------------
    grid_perf <- data.frame(grid, mean_bs = NA_real_, n_ok = 0L)
    
    for (g in 1:nrow(grid)) {
      
      bs_vals <- rep(NA_real_, K_inner)
      ok      <- rep(FALSE,   K_inner)
      
      for (k in 1:K_inner) {
        
        inner_train_idx <- which(folds_inner != k)
        inner_val_idx   <- which(folds_inner == k)
        
        dtr <- dtrain_outer[inner_train_idx, , drop = FALSE]
        dva <- dtrain_outer[inner_val_idx,   , drop = FALSE]
        
        if (sum(dtr[[event_col]] == 1) < 2 || sum(dva[[event_col]] == 1) < 1) next
        
        max_va_time <- max(dva[[time_col]], na.rm = TRUE)
        times_inner <- times[times <= max_va_time]
        if (length(times_inner) == 0) next
        
        # Fit SRF on INNER train (importance off to keep tuning fast)
        set.seed(seed + o*100 + k*10 + g)
        fit_inner <- try(
          rfsrc(
            formula   = srf_formula,
            data      = dtr,
            ntree     = grid$ntree[g],
            mtry      = grid$mtry[g],
            nodesize  = grid$nodesize[g],
            nsplit    = grid$nsplit[g],
            forest    = TRUE,
            importance = "none"
          ),
          silent = TRUE
        )
        if (inherits(fit_inner, "try-error")) next
        
        # Score on INNER validation
        dscore <- dva
        dscore$time  <- dscore[[time_col]]
        dscore$event <- dscore[[event_col]]
        
        score_obj <- try(
          Score(
            object = list("SRF" = fit_inner),
            formula = Surv(time, event) ~ 1,
            data = dscore,
            metrics = "Brier",
            times = times_inner,
            null.model = TRUE,
            summary = "IPA"
          ),
          silent = TRUE
        )
        if (inherits(score_obj, "try-error")) next
        
        bs_model <- score_obj$Brier$score$Brier[
          score_obj$Brier$score$model == "SRF"
        ]
        if (length(bs_model) == 0 || all(!is.finite(bs_model))) next
        
        bs_vals[k] <- mean(bs_model, na.rm = TRUE)  # objective
        ok[k] <- TRUE
      }
      
      grid_perf$mean_bs[g] <- mean(bs_vals, na.rm = TRUE)
      grid_perf$n_ok[g]    <- sum(ok)
    }
    
    inner_grid_tables[[o]] <- grid_perf
    
    ok_rows <- grid_perf[grid_perf$n_ok > 0 & is.finite(grid_perf$mean_bs), , drop = FALSE]
    if (nrow(ok_rows) == 0) {
      if (verbose) message("Outer fold ", o, ": inner tuning failed; skipping.")
      next
    }
    
    best_row <- ok_rows[which.min(ok_rows$mean_bs), , drop = FALSE]
    best_params <- as.list(best_row[, c("mtry", "nodesize", "nsplit", "ntree")])
    best_params_outer[[o]] <- best_params
    
    if (verbose) {
      message("Best params (inner CV): ",
              paste(names(best_params), unlist(best_params), sep = "=", collapse = ", "),
              " | inner mean BS=", round(best_row$mean_bs, 4),
              " | valid inner folds=", best_row$n_ok)
    }
    
    # ----------------------------------------------------------
    # Refit on FULL OUTER train using best params + compute VIMP
    # ----------------------------------------------------------
    set.seed(seed + o)
    fit_outer <- try(
      rfsrc(
        formula   = srf_formula,
        data      = dtrain_outer,
        ntree     = best_params$ntree,
        mtry      = best_params$mtry,
        nodesize  = best_params$nodesize,
        nsplit    = best_params$nsplit,
        forest    = TRUE,
        importance = "permute"   # <-- Option A: fold-wise VIMP
      ),
      silent = TRUE
    )
    if (inherits(fit_outer, "try-error")) {
      if (verbose) message("Outer fold ", o, ": refit failed; skipping.")
      next
    }
    
    # store fold-wise VIMP (named numeric vector)
    vimp_outer[[o]] <- fit_outer$importance
    
    # ----------------------------------------------------------
    # OUTER evaluation on held-out test fold (BS + IPA)
    # ----------------------------------------------------------
    dscore_test <- dtest_outer
    dscore_test$time  <- dscore_test[[time_col]]
    dscore_test$event <- dscore_test[[event_col]]
    
    score_out <- try(
      Score(
        object = list("SRF" = fit_outer),
        formula = Surv(time, event) ~ 1,
        data = dscore_test,
        metrics = "Brier",
        times = times_outer,
        null.model = TRUE,
        summary = "IPA"
      ),
      silent = TRUE
    )
    if (inherits(score_out, "try-error")) {
      if (verbose) message("Outer fold ", o, ": Score failed; skipping.")
      next
    }
    
    bs_model <- score_out$Brier$score$Brier[score_out$Brier$score$model == "SRF"]
    bs_null  <- score_out$Brier$score$Brier[score_out$Brier$score$model == "Null model"]
    
    bs_outer[o]      <- mean(bs_model, na.rm = TRUE)
    bs_null_outer[o] <- mean(bs_null,  na.rm = TRUE)
    ipa_outer[o]     <- 1 - (bs_outer[o] / bs_null_outer[o])
    
    if (verbose) {
      message("OUTER fold ", o,
              " | BS=", round(bs_outer[o], 4),
              " | BS_null=", round(bs_null_outer[o], 4),
              " | IPA=", round(ipa_outer[o], 4))
    }
  }
  
  # -----------------------------
  # Summarize across outer folds
  # -----------------------------
  summary_df <- data.frame(
    rows_used = nrow(dsub),
    events = sum(dsub[[event_col]] == 1),
    folds_scored = sum(is.finite(bs_outer)),
    BS_mean = mean(bs_outer, na.rm = TRUE),
    BS_sd   = sd(bs_outer,   na.rm = TRUE),
    BS_null_mean = mean(bs_null_outer, na.rm = TRUE),
    BS_null_sd   = sd(bs_null_outer,   na.rm = TRUE),
    IPA_mean = mean(ipa_outer, na.rm = TRUE),
    IPA_sd   = sd(ipa_outer,   na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  list(
    summary = summary_df,
    fold_results = data.frame(
      fold = 1:K_outer,
      BS = bs_outer,
      BS_null = bs_null_outer,
      IPA = ipa_outer
    ),
    best_params_outer = best_params_outer,
    inner_grid_tables = inner_grid_tables,
    vimp_outer = vimp_outer
  )
}


#

run_doublecv_srf_over_combos <- function(
    df,
    var_groups,
    time_col    = "tfs_months",
    event_col   = "event_death",
    combo_sizes = 1:length(var_groups),
    
    # pass-through to the working RSF double-CV
    K_outer = 10,
    K_inner = 5,
    times   = seq(12, 60, 12),
    grid    = expand.grid(
      mtry     = c(2, 4, 8),
      nodesize = c(5, 10, 15),
      nsplit   = c(10),
      ntree    = c(500),
      stringsAsFactors = FALSE
    ),
    
    # optional filters like lasso pipeline
    min_rows   = 30,
    min_events = 10,
    
    seed    = 123,
    verbose = TRUE
) {
  
  # generate combos of group names (e.g., "lab", "lab+ba", ...)
  combos <- make_var_group_combos(var_groups, sizes = combo_sizes)
  
  summary_rows <- list()
  details_rows <- list()
  
  combo_id <- 1L
  
  for (groups_in in combos) {
    
    combo_label <- paste(groups_in, collapse = "+")
    vars <- unique(unlist(var_groups[groups_in], use.names = FALSE))
    
    if (verbose) {
      message("\n==============================")
      message("COMBO ", combo_id, ": ", combo_label, " (n_vars=", length(vars), ")")
      message("==============================")
    }
    
    # check columns exist
    needed <- c(time_col, event_col, vars)
    missing_cols <- setdiff(needed, colnames(df))
    if (length(missing_cols) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        combo_id = combo_id,
        groups   = combo_label,
        n_vars   = length(vars),
        rows_used = NA_integer_,
        events    = NA_integer_,
        folds_scored = 0L,
        BS_mean = NA_real_, BS_sd = NA_real_,
        BS_null_mean = NA_real_, BS_null_sd = NA_real_,
        IPA_mean = NA_real_, IPA_sd = NA_real_,
        status = paste("missing:", paste(missing_cols, collapse = ",")),
        stringsAsFactors = FALSE
      )
      details_rows[[as.character(combo_id)]] <- list(
        groups = groups_in,
        vars   = vars,
        result = NULL
      )
      combo_id <- combo_id + 1L
      next
    }
    
    # quick sufficiency check on complete cases (combo-specific)
    dsub <- df[complete.cases(df[, needed]), needed, drop = FALSE]
    n_rows <- nrow(dsub)
    n_events <- sum(dsub[[event_col]] == 1)
    
    if (n_rows < min_rows || n_events < min_events) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        combo_id = combo_id,
        groups   = combo_label,
        n_vars   = length(vars),
        rows_used = n_rows,
        events    = n_events,
        folds_scored = 0L,
        BS_mean = NA_real_, BS_sd = NA_real_,
        BS_null_mean = NA_real_, BS_null_sd = NA_real_,
        IPA_mean = NA_real_, IPA_sd = NA_real_,
        status = "insufficient rows/events",
        stringsAsFactors = FALSE
      )
      details_rows[[as.character(combo_id)]] <- list(
        groups = groups_in,
        vars   = vars,
        result = NULL
      )
      combo_id <- combo_id + 1L
      next
    }
    
    # ---- RUN THE WORKING FUNCTION (engine) ----
    res <- doublecv_srf_brier(
      df        = df,
      vars      = vars,
      time_col  = time_col,
      event_col = event_col,
      K_outer   = K_outer,
      K_inner   = K_inner,
      times     = times,
      grid      = grid,
      seed      = seed,
      verbose   = verbose
    )
    
    # collect summary row (like lasso summary table)
    s <- res$summary
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      combo_id = combo_id,
      groups   = combo_label,
      n_vars   = length(vars),
      rows_used = s$rows_used,
      events    = s$events,
      folds_scored = s$folds_scored,
      BS_mean = s$BS_mean,
      BS_sd   = s$BS_sd,
      BS_null_mean = s$BS_null_mean,
      BS_null_sd   = s$BS_null_sd,
      IPA_mean = s$IPA_mean,
      IPA_sd   = s$IPA_sd,
      status   = "ok",
      stringsAsFactors = FALSE
    )
    
    details_rows[[as.character(combo_id)]] <- list(
      groups = groups_in,
      vars   = vars,
      result = res
    )
    
    combo_id <- combo_id + 1L
  }
  
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  
  list(
    summary = summary_df,
    details = details_rows
  )
}

var_groups <- list(
  lab      = lab_vars,
  ba       = ba_vars,
  ram      = ram_vars,
  clinical = clinical_vars
)

res_srf_all <- run_doublecv_srf_over_combos(
  df = df_imputed,
  var_groups = var_groups,
  combo_sizes = 1:4,
  K_outer = 10,
  K_inner = 5,
  times = seq(12, 60, 12),
  grid = expand.grid(
    mtry = c(2,4,8),
    nodesize = c(5,10,15),
    nsplit = c(10),
    ntree = c(500),
    stringsAsFactors = FALSE
  ),
  min_rows = 150,
  min_events = 10,
  seed = 123,
  verbose = TRUE
)

res_srf_all$summary





#-----------------------------------------------------
#Ensemble 2 (Tree--Cluster Cox Model)
#-----------------------------------------------------



#A helper function. Terminal-node assignment for rpart using type="matrix"
#Given a fitted rpart tree and new data, return the terminal node ID
#   (a.k.a. the leaf node) for each new observation.

.rpart_terminal_from_predict_matrix <- function(tree_fit, newdata) {
  # predict(..., type="matrix") returns a matrix where each row corresponds
  # to an observation, and columns correspond to node IDs along the path.
  # Values become NA after the observation terminates (reaches a leaf).
  pm <- predict(tree_fit, newdata = newdata, type = "matrix")
  apply(pm, 1, function(z) {
    z <- z[!is.na(z)]
    z[length(z)]
  })
}

# ------------------------------------------------------------
# Stage 2: build tree + Cox-per-cluster given top_vars
# ------------------------------------------------------------
build_ens2_stage2_tree_clustercox <- function(
    df,
    top_vars,
    time_col  = "tfs_months",
    event_col = "event_death",
    #Tuned in outerRank function
    tree_maxdepth  = 2,
    #Fixed hyperparams
    tree_minsplit  = 30,
    tree_minbucket = 15,
    tree_cp        = 0.01,
    #   if a cluster has fewer than this many events => skip Cox for that cluster
    min_cluster_events = 5,
    
    seed = 123,
    verbose = FALSE
) {
  
  #clean + check,  Ensure top_vars are unique, trimmed, non-empty, and exist in df.
  top_vars <- unique(trimws(top_vars))
  top_vars <- top_vars[top_vars != ""]
  top_vars <- top_vars[top_vars %in% colnames(df)]
  if (length(top_vars) == 0) stop("No valid top_vars found in df.")
  
  # complete-case subset for stage2 
  needed <- c(time_col, event_col, top_vars)
  dsub <- df[stats::complete.cases(df[, needed, drop = FALSE]), needed, drop = FALSE]
  if (nrow(dsub) == 0) stop("No rows left after complete.cases().")
  dsub[[time_col]] <- ifelse(dsub[[time_col]] <= 0, 1e-6, dsub[[time_col]])
  
  safe_top <- paste0("`", top_vars, "`")
  
  # ---- tree for clustering ----
  f_tree <- stats::as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(safe_top, collapse = " + "))
  )
  
  # Fit tree: model=TRUE: keep model.frame inside tree_fit (helps prediction stability)
  #   - control: depth and minimum split sizes control number of clusters/leaves
  set.seed(seed)
  tree_fit <- rpart::rpart(
    formula = f_tree,
    data    = dsub,
    method  = "exp",
    model   = TRUE,
    control = rpart::rpart.control(
      maxdepth  = tree_maxdepth,
      minsplit  = tree_minsplit,
      minbucket = tree_minbucket,
      cp        = tree_cp,
      xval      = 0
    )
  )
  
  # assign clusters for training
  dsub$cluster <- factor(tree_fit$where)
  
  # Cluster diagnostics table:
  #   - n: number of samples per cluster
  #   - events: number of events per cluster
  cluster_tab <- data.frame(
    cluster = levels(dsub$cluster),
    n      = as.integer(table(dsub$cluster)),
    events = as.integer(tapply(dsub[[event_col]], dsub$cluster, sum)),
    stringsAsFactors = FALSE
  )
  
  
  #Now we fit cox model on each clusters on the same top vars 
  # Cox per cluster
  f_cox <- stats::as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(safe_top, collapse = " + "))
  )
  
  cox_models <- list()
  dropped_clusters <- character(0)
  
  for (cl in levels(dsub$cluster)) {
    dcl <- dsub[dsub$cluster == cl, , drop = FALSE]
    n_ev <- sum(dcl[[event_col]] == 1)
    
    if (n_ev < min_cluster_events) {
      dropped_clusters <- c(dropped_clusters, cl)
      next
    }
    #fit cox 
    fit_cl <- try(survival::coxph(f_cox, data = dcl, x = TRUE), silent = TRUE)
    if (inherits(fit_cl, "try-error")) {
      dropped_clusters <- c(dropped_clusters, cl)
      next
    }
    cox_models[[as.character(cl)]] <- fit_cl
  }
  
  # global fallback Cox (fit once per stage2 model), only if the test data lands in cluster with n cluster specific cox model
  cox_global <- try(survival::coxph(f_cox, data = dsub, x = TRUE), silent = TRUE)
  
  #  prediction function , tree cluster -> corresponding Cox -> survival
  predict_stage2 <- function(newdata, times) {
    miss <- setdiff(top_vars, colnames(newdata))
    if (length(miss) > 0) stop("newdata missing predictors: ", paste(miss, collapse = ","))
    
    times <- sort(unique(times))
    # no. of test obs
    n_new <- nrow(newdata)
    
    # terminal node IDs robustly
    # Use the helper to get terminal-node IDs for each new row.
    cl_new <- .rpart_terminal_from_predict_matrix(tree_fit, newdata)
    cl_new <- as.character(cl_new)
    
    # Allocate survival matrix:
    #   rows: new observations
    #   cols: requested times
    S_mat <- matrix(NA_real_, nrow = n_new, ncol = length(times))
    colnames(S_mat) <- paste0("t=", times)
    
    #Another helper to fill survival rows from a iven cox model
    fill_from_cox <- function(fit, idx) {
      if (length(idx) == 0) return()
      # baseline surv S_0(t) from cox fit
      # survfit(fit) without newdata gives baseline.
      base <- survival::survfit(fit)
      
      # Turn baseline survival into a step function:
      S0_at <- stats::stepfun(base$time, c(1, base$surv))
      # Linear predictor lp = x^T beta for each row, using Cox fit
      lp <- try(stats::predict(fit,
                               newdata = newdata[idx, top_vars, drop = FALSE],
                               type = "lp"),
                silent = TRUE)
      if (inherits(lp, "try-error")) return()
      
      # Cox survival formula:
      #   S(t|x) = S0(t) ^ exp(lp)
      exp_lp <- exp(lp)
      for (j in seq_along(idx)) {
        S_mat[idx[j], ] <<- S0_at(times) ^ exp_lp[j]
      }
    }
    # Fill survival matrix cluster-by-cluster
    for (cl in unique(cl_new)) {
      idx <- which(cl_new == cl)
      fit_cl <- cox_models[[cl]]
      
      if (!is.null(fit_cl)) {
        fill_from_cox(fit_cl, idx)
      } else if (!inherits(cox_global, "try-error")) {
        fill_from_cox(cox_global, idx)
      }
    }
    
    # final safety: fill any NA with KM survival from training
    # riskRegression::Score cannot handle NA risk/survival.
    # If anything is NA due to model failures, fill with KM survival.
    if (anyNA(S_mat)) {
      km <- survival::survfit(survival::Surv(dsub[[time_col]], dsub[[event_col]]) ~ 1)
      S_km <- stats::stepfun(km$time, c(1, km$surv))(times)
      # rows with any non-finite values
      bad <- which(apply(S_mat, 1, function(z) any(!is.finite(z))))
      # Fill those rows with KM survival
      for (i in bad) S_mat[i, ] <- S_km
    }
    # Return survival + risk matrix.
    list(cluster = cl_new, survival = S_mat, risk = 1 - S_mat)
  }
  
  list(
    data_info = list(rows_used = nrow(dsub), events = sum(dsub[[event_col]] == 1)),
    top_vars  = top_vars,
    tree = list(fit = tree_fit, cluster_table = cluster_tab),
    cox_by_cluster = list(models = cox_models, dropped_clusters = dropped_clusters, formula = f_cox),
    predict = predict_stage2
  )
}

#double cv Function: outer RSF ranking once; inner tunes stage2 only
doublecv_ensemble2_brier_outerRank <- function(
    df,
    vars,
    time_col  = "tfs_months",
    event_col = "event_death",
    K_outer   = 10,
    K_inner   = 5,
    times     = seq(12, 60, 12),
    
    grid = expand.grid(
      top_n = c(5, 10, 15),
      tree_maxdepth = c(1, 2, 3),
      stringsAsFactors = FALSE
    ),
    
    # RSF used ONLY for outer ranking
    ntree_rank = 500,
    mtry      = NULL,
    nodesize  = 15,
    nsplit    = 10,
    importance_rank = "permute",
    
    # stage2 fixed params
    tree_minsplit  = 30,
    tree_minbucket = 15,
    tree_cp        = 0.01,
    min_cluster_events = 5,
    
    seed    = 123,
    verbose = TRUE
) {
  if (!requireNamespace("riskRegression", quietly = TRUE)) stop("Install riskRegression (Score)")
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) stop("Install randomForestSRC")
  if (!requireNamespace("survival", quietly = TRUE)) stop("Install survival")
  suppressPackageStartupMessages(library(survival))
  
  # complete-case subset on ALL candidate vars (ranking needs them)
  needed <- c(time_col, event_col, vars)
  dsub <- df[stats::complete.cases(df[, needed, drop = FALSE]), needed, drop = FALSE]
  if (nrow(dsub) == 0) stop("No rows left after complete.cases().")
  dsub[[time_col]] <- ifelse(dsub[[time_col]] <= 0, 1e-6, dsub[[time_col]])
  dsub <- dsub[order(rownames(dsub)), , drop = FALSE]
  
  # outer folds (stratified)
  set.seed(seed)
  n <- nrow(dsub)
  e <- as.integer(dsub[[event_col]])
  folds_outer <- integer(n)
  idx1 <- which(e == 1); idx0 <- which(e == 0)
  folds_outer[idx1] <- sample(rep(1:K_outer, length.out = length(idx1)))
  folds_outer[idx0] <- sample(rep(1:K_outer, length.out = length(idx0)))
  
  # storage
  bs_outer <- bs_null_outer <- ipa_outer <- rep(NA_real_, K_outer)
  
  best_params_outer <- vector("list", K_outer)
  inner_grid_tables <- vector("list", K_outer)
  
  ranked_vars_outer <- vector("list", K_outer)   # full ranking per outer fold
  top_vars_outer    <- vector("list", K_outer)   # chosen top_vars per outer fold
  rank_vimp_outer   <- vector("list", K_outer)   # (optional) VIMP vector
  
  cluster_tables_outer   <- vector("list", K_outer)
  dropped_clusters_outer <- vector("list", K_outer)
  
  surv_formula_score <- Surv(time, event) ~ 1
  
  # default mtry
  if (is.null(mtry)) mtry <- max(1, floor(sqrt(length(vars))))
  safe_vars <- paste0("`", vars, "`")
  
  # RSF formula
  f_rsf <- stats::as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(safe_vars, collapse = " + "))
  )
  
  for (o in 1:K_outer) {
    if (verbose) message("\n===== OUTER FOLD ", o, "/", K_outer, " (Outer-Rank Ensemble2) =====")
    
    train_idx <- which(folds_outer != o)
    test_idx  <- which(folds_outer == o)
    
    dtrain_outer <- dsub[train_idx, , drop = FALSE]
    dtest_outer  <- dsub[test_idx,  , drop = FALSE]
    
    if (sum(dtrain_outer[[event_col]] == 1) < 2 || sum(dtest_outer[[event_col]] == 1) < 1) {
      if (verbose) message("Outer fold ", o, ": too few events; skipping.")
      next
    }
    
    max_test_time <- max(dtest_outer[[time_col]], na.rm = TRUE)
    times_outer <- times[times <= max_test_time]
    if (length(times_outer) == 0) {
      if (verbose) message("Outer fold ", o, ": no usable eval times; skipping.")
      next
    }
    
    # ----------------------------------------------------------
    # 1) Fit ONE RSF on OUTER train for ranking of VIMP
    # ----------------------------------------------------------
    set.seed(seed + o)
    fit_rank <- try(
      randomForestSRC::rfsrc(
        formula    = f_rsf,
        data       = dtrain_outer,
        ntree      = ntree_rank,
        mtry       = mtry,
        nodesize   = nodesize,
        nsplit     = nsplit,
        forest     = TRUE,
        importance = importance_rank
      ),
      silent = TRUE
    )
    if (inherits(fit_rank, "try-error")) {
      if (verbose) message("Outer fold ", o, ": ranking RSF failed; skipping.")
      next
    }
    
    vimp <- fit_rank$importance
    if (is.null(vimp) || length(vimp) == 0) {
      if (verbose) message("Outer fold ", o, ": ranking RSF returned no importance; skipping.")
      next
    }
    vimp_sorted <- sort(vimp, decreasing = TRUE)
    ranked_vars <- names(vimp_sorted)
    
    ranked_vars_outer[[o]] <- ranked_vars
    rank_vimp_outer[[o]]   <- vimp_sorted
    
    # ----------------------------------------------------------
    # 2) INNER CV: tune using ONLY stage2 fits (no RSF)
    # ----------------------------------------------------------
    set.seed(seed + o)
    ntr <- nrow(dtrain_outer)
    e_tr <- as.integer(dtrain_outer[[event_col]])
    folds_inner <- integer(ntr)
    id1 <- which(e_tr == 1); id0 <- which(e_tr == 0)
    folds_inner[id1] <- sample(rep(1:K_inner, length.out = length(id1)))
    folds_inner[id0] <- sample(rep(1:K_inner, length.out = length(id0)))
    
    # Heres the INNER tuning over grid(top_n, tree_maxdepth). minimize mean Brier score
    grid_perf <- data.frame(grid, mean_bs = NA_real_, n_ok = 0L)
    
    for (g in 1:nrow(grid)) {
      bs_vals <- rep(NA_real_, K_inner)
      ok <- rep(FALSE, K_inner)
      
      top_n_g <- grid$top_n[g]
      depth_g <- grid$tree_maxdepth[g]
      top_vars_g <- ranked_vars[seq_len(min(top_n_g, length(ranked_vars)))]
      
      for (k in 1:K_inner) {
        inner_train_idx <- which(folds_inner != k)
        inner_val_idx   <- which(folds_inner == k)
        
        dtr <- dtrain_outer[inner_train_idx, , drop = FALSE]
        dva <- dtrain_outer[inner_val_idx,   , drop = FALSE]
        
        if (sum(dtr[[event_col]] == 1) < 2 || sum(dva[[event_col]] == 1) < 1) next
        
        max_va_time <- max(dva[[time_col]], na.rm = TRUE)
        times_inner <- times[times <= max_va_time]
        if (length(times_inner) == 0) next
        
        # stage2 fit ONLY (tree+cox) using candidate top_vars + depth
        st2 <- try(
          build_ens2_stage2_tree_clustercox(
            df = dtr,
            top_vars = top_vars_g,
            time_col = time_col,
            event_col = event_col,
            tree_maxdepth = depth_g,
            tree_minsplit = tree_minsplit,
            tree_minbucket = tree_minbucket,
            tree_cp = tree_cp,
            min_cluster_events = min_cluster_events,
            seed = seed + o*100 + g*10 + k,
            verbose = FALSE
          ),
          silent = TRUE
        )
        if (inherits(st2, "try-error")) next
        
        pred <- try(st2$predict(dva, times = times_inner), silent = TRUE)
        if (inherits(pred, "try-error")) next
        
        risk_mat <- pred$risk
        colnames(risk_mat) <- times_inner
        
        dscore <- dva
        dscore$time  <- dscore[[time_col]]
        dscore$event <- dscore[[event_col]]
        
        #compute the BS score
        sc <- try(
          riskRegression::Score(
            object = list("Ensemble2" = risk_mat),
            formula = Surv(time, event) ~ 1,
            data = dscore,
            metrics = "Brier",
            times = times_inner,
            null.model = TRUE,
            summary = "IPA"
          ),
          silent = TRUE
        )
        if (inherits(sc, "try-error")) next
        
        bs_model <- sc$Brier$score$Brier[sc$Brier$score$model == "Ensemble2"]
        if (length(bs_model) == 0 || all(!is.finite(bs_model))) next
        
        bs_vals[k] <- mean(bs_model, na.rm = TRUE)
        ok[k] <- TRUE
      }
      #Store results
      grid_perf$mean_bs[g] <- mean(bs_vals, na.rm = TRUE)
      grid_perf$n_ok[g]    <- sum(ok)
    }
    
    inner_grid_tables[[o]] <- grid_perf
    
    ok_rows <- grid_perf[grid_perf$n_ok > 0 & is.finite(grid_perf$mean_bs), , drop = FALSE]
    if (nrow(ok_rows) == 0) {
      if (verbose) message("Outer fold ", o, ": inner tuning failed; skipping.")
      next
    }
    
    best_row <- ok_rows[which.min(ok_rows$mean_bs), , drop = FALSE]
    best_params <- as.list(best_row[, c("top_n", "tree_maxdepth")])
    best_params_outer[[o]] <- best_params
    
    # chosen top_vars on outer fold
    top_vars_best <- ranked_vars[seq_len(min(best_params$top_n, length(ranked_vars)))]
    top_vars_outer[[o]] <- top_vars_best
    
    if (verbose) {
      message("Best params: top_n=", best_params$top_n,
              ", tree_maxdepth=", best_params$tree_maxdepth,
              " | inner mean BS=", round(best_row$mean_bs, 4),
              " | valid inner folds=", best_row$n_ok)
    }
    
    # ----------------------------------------------------------
    # 3) Refit stage2 on FULL OUTER train with best params
    # ----------------------------------------------------------
    st2_outer <- try(
      build_ens2_stage2_tree_clustercox(
        df = dtrain_outer,
        top_vars = top_vars_best,
        time_col = time_col,
        event_col = event_col,
        tree_maxdepth = best_params$tree_maxdepth,
        tree_minsplit = tree_minsplit,
        tree_minbucket = tree_minbucket,
        tree_cp = tree_cp,
        min_cluster_events = min_cluster_events,
        seed = seed + o,
        verbose = FALSE
      ),
      silent = TRUE
    )
    if (inherits(st2_outer, "try-error")) {
      if (verbose) message("Outer fold ", o, ": stage2 refit failed; skipping.")
      next
    }
    
    cluster_tables_outer[[o]]   <- st2_outer$tree$cluster_table
    dropped_clusters_outer[[o]] <- st2_outer$cox_by_cluster$dropped_clusters
    
    # ----------------------------------------------------------
    # 4) OUTER evaluation
    # ----------------------------------------------------------
    pred_out <- try(st2_outer$predict(dtest_outer, times = times_outer), silent = TRUE)
    if (inherits(pred_out, "try-error")) {
      if (verbose) message("Outer fold ", o, ": predict failed; skipping.")
      next
    }
    
    risk_mat <- pred_out$risk
    colnames(risk_mat) <- times_outer
    
    dscore_test <- dtest_outer
    dscore_test$time  <- dscore_test[[time_col]]
    dscore_test$event <- dscore_test[[event_col]]
    
    sc_out <- try(
      riskRegression::Score(
        object = list("Ensemble2" = risk_mat),
        formula = surv_formula_score,
        data = dscore_test,
        metrics = "Brier",
        times = times_outer,
        null.model = TRUE,
        summary = "IPA"
      ),
      silent = TRUE
    )
    if (inherits(sc_out, "try-error")) {
      if (verbose) message("Outer fold ", o, ": Score failed; skipping.")
      next
    }
    # Extract Brier curves for model and null model
    bs_model <- sc_out$Brier$score$Brier[sc_out$Brier$score$model == "Ensemble2"]
    bs_null  <- sc_out$Brier$score$Brier[sc_out$Brier$score$model == "Null model"]
    
    bs_outer[o]      <- mean(bs_model, na.rm = TRUE)
    bs_null_outer[o] <- mean(bs_null,  na.rm = TRUE)
    ipa_outer[o]     <- 1 - (bs_outer[o] / bs_null_outer[o])
    
    if (verbose) {
      message("OUTER fold ", o,
              " | BS=", round(bs_outer[o], 4),
              " | BS_null=", round(bs_null_outer[o], 4),
              " | IPA=", round(ipa_outer[o], 4))
    }
  }
  #Summarize across outer folds
  summary_df <- data.frame(
    rows_used = nrow(dsub),
    events = sum(dsub[[event_col]] == 1),
    folds_scored = sum(is.finite(bs_outer)),
    BS_mean = mean(bs_outer, na.rm = TRUE),
    BS_sd   = sd(bs_outer,   na.rm = TRUE),
    BS_null_mean = mean(bs_null_outer, na.rm = TRUE),
    BS_null_sd   = sd(bs_null_outer,   na.rm = TRUE),
    IPA_mean = mean(ipa_outer, na.rm = TRUE),
    IPA_sd   = sd(ipa_outer,   na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  
  fold_results <- data.frame(
    fold = 1:K_outer,
    BS = bs_outer,
    BS_null = bs_null_outer,
    IPA = ipa_outer,
    top_n = sapply(best_params_outer, function(x) if (is.null(x)) NA else x$top_n),
    tree_maxdepth = sapply(best_params_outer, function(x) if (is.null(x)) NA else x$tree_maxdepth),
    stringsAsFactors = FALSE
  )
  
  list(
    summary = summary_df,
    fold_results = fold_results,
    best_params_outer = best_params_outer,
    inner_grid_tables = inner_grid_tables,
    
    ranked_vars_outer = ranked_vars_outer,
    rank_vimp_outer   = rank_vimp_outer,
    top_vars_outer    = top_vars_outer,
    
    cluster_tables_outer = cluster_tables_outer,
    dropped_clusters_outer = dropped_clusters_outer
  )
}




# ------------------------------------------------------------
# Helpers: combos generator (reuse the existing one)
# ------------------------------------------------------------
make_var_group_combos <- function(var_groups, sizes = 1:length(var_groups)) {
  group_names <- names(var_groups)
  combos <- list()
  id <- 1
  for (k in sizes) {
    cmb <- utils::combn(group_names, k, simplify = FALSE)
    for (g in cmb) {
      combos[[id]] <- g
      id <- id + 1
    }
  }
  combos
}

# ------------------------------------------------------------
# Runner: double-CV Ensemble2 over all var-group combos
# Stores tuned params + top_vars per outer fold via the engine
run_doublecv_ens2_over_combos <- function(
    df,
    var_groups,
    time_col    = "tfs_months",
    event_col   = "event_death",
    combo_sizes = 1:length(var_groups),
    
    # pass-through to engine
    K_outer = 10,
    K_inner = 5,
    times   = seq(12, 60, 12),
    
    grid = expand.grid(
      top_n = c(5, 10, 15),
      tree_maxdepth = c(1, 2, 3),
      stringsAsFactors = FALSE
    ),
    
    # RSF ranking params (once per outer fold)
    ntree_rank = 500,
    mtry      = NULL,
    nodesize  = 15,
    nsplit    = 10,
    importance_rank = "anti",
    
    # stage2 fixed params
    tree_minsplit  = 30,
    tree_minbucket = 15,
    tree_cp        = 0.01,
    min_cluster_events = 5,
    
    # filters
    min_rows   = 30,
    min_events = 10,
    
    seed    = 123,
    verbose = TRUE
) {
  combos <- make_var_group_combos(var_groups, sizes = combo_sizes)
  
  summary_rows <- list()
  details_rows <- list()
  combo_id <- 1L
  
  for (groups_in in combos) {
    combo_label <- paste(groups_in, collapse = "+")
    vars <- unique(unlist(var_groups[groups_in], use.names = FALSE))
    
    if (verbose) {
      message("\n==============================")
      message("COMBO ", combo_id, ": ", combo_label, " (n_vars=", length(vars), ")")
      message("==============================")
    }
    
    needed <- c(time_col, event_col, vars)
    missing_cols <- setdiff(needed, colnames(df))
    if (length(missing_cols) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        combo_id = combo_id,
        groups   = combo_label,
        n_vars   = length(vars),
        rows_used = NA_integer_,
        events    = NA_integer_,
        folds_scored = 0L,
        BS_mean = NA_real_, BS_sd = NA_real_,
        BS_null_mean = NA_real_, BS_null_sd = NA_real_,
        IPA_mean = NA_real_, IPA_sd = NA_real_,
        status = paste("missing:", paste(missing_cols, collapse = ",")),
        stringsAsFactors = FALSE
      )
      details_rows[[as.character(combo_id)]] <- list(groups = groups_in, vars = vars, result = NULL)
      combo_id <- combo_id + 1L
      next
    }
    
    dcheck <- df[stats::complete.cases(df[, needed, drop = FALSE]), needed, drop = FALSE]
    n_rows   <- nrow(dcheck)
    n_events <- sum(dcheck[[event_col]] == 1)
    
    if (n_rows < min_rows || n_events < min_events) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        combo_id = combo_id,
        groups   = combo_label,
        n_vars   = length(vars),
        rows_used = n_rows,
        events    = n_events,
        folds_scored = 0L,
        BS_mean = NA_real_, BS_sd = NA_real_,
        BS_null_mean = NA_real_, BS_null_sd = NA_real_,
        IPA_mean = NA_real_, IPA_sd = NA_real_,
        status = "insufficient rows/events",
        stringsAsFactors = FALSE
      )
      details_rows[[as.character(combo_id)]] <- list(groups = groups_in, vars = vars, result = NULL)
      combo_id <- combo_id + 1L
      next
    }
    
    #run the outer Rank double cv function
    res <- doublecv_ensemble2_brier_outerRank(
      df = df,
      vars = vars,
      time_col = time_col,
      event_col = event_col,
      K_outer = K_outer,
      K_inner = K_inner,
      times = times,
      grid = grid,
      
      ntree_rank = ntree_rank,
      mtry = mtry,
      nodesize = nodesize,
      nsplit = nsplit,
      importance_rank = importance_rank,
      
      tree_minsplit = tree_minsplit,
      tree_minbucket = tree_minbucket,
      tree_cp = tree_cp,
      min_cluster_events = min_cluster_events,
      
      seed = seed,
      verbose = verbose
    )
    
    s <- res$summary
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      combo_id = combo_id,
      groups   = combo_label,
      n_vars   = length(vars),
      rows_used = s$rows_used,
      events    = s$events,
      folds_scored = s$folds_scored,
      BS_mean = s$BS_mean,
      BS_sd   = s$BS_sd,
      BS_null_mean = s$BS_null_mean,
      BS_null_sd   = s$BS_null_sd,
      IPA_mean = s$IPA_mean,
      IPA_sd   = s$IPA_sd,
      status   = "ok",
      stringsAsFactors = FALSE
    )
    
    details_rows[[as.character(combo_id)]] <- list(groups = groups_in, vars = vars, result = res)
    combo_id <- combo_id + 1L
  }
  
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  
  list(summary = summary_df, details = details_rows)
}

var_groups <- list(
  lab      = lab_vars,
  ba       = ba_vars,
  ram      = ram_vars,
  clinical = clinical_vars
)

res_ens2_all <- run_doublecv_ens2_over_combos(
  df = df_imputed,
  var_groups = var_groups,
  combo_sizes = 1:4,
  K_outer = 10,
  K_inner = 5,
  times = seq(12, 60, 12),
  grid = expand.grid(
    top_n = c(5, 10, 15),
    tree_maxdepth = c(1, 2),
    stringsAsFactors = FALSE
  ),
  min_rows = 150,
  min_events = 10,
  seed = 123,
  verbose = TRUE
)
res_ens2_all$summary

#HYbrid

doublecv_hybrid_srf_topn_cox_brier <- function(
    df,
    vars,
    time_col  = "tfs_months",
    event_col = "event_death",
    K_outer   = 10,
    K_inner   = 5,
    times     = seq(12, 60, 12),
    
    # tune only this
    grid_top_n = c(5, 10, 15, 20),
    
    # SRF ranking (outer only)
    ntree_rank = 500,
    mtry       = NULL,
    nodesize   = 15,
    nsplit     = 10,
    importance_rank = "anti",
    
    seed    = 123,
    verbose = TRUE
) {
  if (!requireNamespace("riskRegression", quietly = TRUE)) stop("Install riskRegression")
  if (!requireNamespace("randomForestSRC", quietly = TRUE)) stop("Install randomForestSRC")
  if (!requireNamespace("survival", quietly = TRUE)) stop("Install survival")
  suppressPackageStartupMessages(library(survival))
  
  # ----------------------------------------------------------
  # INTERNAL helper: stratified folds
  # ----------------------------------------------------------
  make_stratified_folds <- function(event, K, seed) {
    set.seed(seed)
    e <- as.integer(event)
    n <- length(e)
    folds <- integer(n)
    
    idx1 <- which(e == 1)
    idx0 <- which(e == 0)
    
    folds[idx1] <- sample(rep(1:K, length.out = length(idx1)))
    folds[idx0] <- sample(rep(1:K, length.out = length(idx0)))
    
    folds
  }
  
  # ----------------------------------------------------------
  # clean vars
  # ----------------------------------------------------------
  vars <- unique(trimws(vars))
  vars <- vars[vars != ""]
  vars <- vars[vars %in% colnames(df)]
  if (length(vars) == 0) stop("No valid predictors.")
  
  needed <- c(time_col, event_col, vars)
  dsub <- df[stats::complete.cases(df[, needed, drop = FALSE]), needed, drop = FALSE]
  if (nrow(dsub) == 0) stop("No rows after complete.cases().")
  dsub[[time_col]] <- ifelse(dsub[[time_col]] <= 0, 1e-6, dsub[[time_col]])
  dsub <- dsub[order(rownames(dsub)), , drop = FALSE]
  
  if (is.null(mtry)) mtry <- max(1, floor(sqrt(length(vars))))
  
  safe_vars <- paste0("`", vars, "`")
  f_rsf <- stats::as.formula(
    paste0("Surv(", time_col, ", ", event_col, ") ~ ", paste(safe_vars, collapse = " + "))
  )
  
  # ----------------------------------------------------------
  # OUTER folds
  # ----------------------------------------------------------
  folds_outer <- make_stratified_folds(dsub[[event_col]], K_outer, seed)
  
  bs_outer      <- rep(NA_real_, K_outer)
  bs_null_outer <- rep(NA_real_, K_outer)
  ipa_outer     <- rep(NA_real_, K_outer)
  
  best_topn_outer   <- rep(NA_integer_, K_outer)
  inner_grid_tables <- vector("list", K_outer)
  
  ranked_vars_outer <- vector("list", K_outer)
  vimp_outer        <- vector("list", K_outer)
  top_vars_outer    <- vector("list", K_outer)
  
  surv_formula_score <- Surv(time, event) ~ 1
  
  # ==========================================================
  # OUTER LOOP
  # ==========================================================
  for (o in 1:K_outer) {
    if (verbose) message("\n===== OUTER FOLD ", o, "/", K_outer, " (Hybrid SRF→Cox) =====")
    
    train_idx <- which(folds_outer != o)
    test_idx  <- which(folds_outer == o)
    
    dtrain_outer <- dsub[train_idx, , drop = FALSE]
    dtest_outer  <- dsub[test_idx,  , drop = FALSE]
    
    n_tr <- nrow(dtrain_outer)
    n_te <- nrow(dtest_outer)
    ev_tr <- sum(dtrain_outer[[event_col]] == 1, na.rm = TRUE)
    ev_te <- sum(dtest_outer[[event_col]] == 1, na.rm = TRUE)
    
    if (verbose) {
      message("Outer fold ", o, ": n_train=", n_tr, " (events=", ev_tr, ")",
              " | n_test=", n_te, " (events=", ev_te, ")")
    }
    
    if (ev_tr < 2 || ev_te < 1) {
      if (verbose) message("Outer fold ", o, ": too few events; skipping.")
      next
    }
    
    max_test_time <- max(dtest_outer[[time_col]], na.rm = TRUE)
    times_outer <- times[times <= max_test_time]
    if (length(times_outer) == 0) {
      if (verbose) message("Outer fold ", o, ": no usable evaluation times (max_test_time=", round(max_test_time, 3), "); skipping.")
      next
    }
    if (verbose) {
      message("Outer fold ", o, ": eval times = [", paste(times_outer, collapse = ", "), "]")
    }
    
    # -----------------------------
    # (1) SRF ranking (outer only)
    # -----------------------------
    if (verbose) message("Outer fold ", o, ": fitting SRF for ranking (importance=", importance_rank, ") …")
    
    set.seed(seed + o)
    fit_rank <- try(
      randomForestSRC::rfsrc(
        formula    = f_rsf,
        data       = dtrain_outer,
        ntree      = ntree_rank,
        mtry       = mtry,
        nodesize   = nodesize,
        nsplit     = nsplit,
        forest     = TRUE,
        importance = importance_rank
      ),
      silent = TRUE
    )
    if (inherits(fit_rank, "try-error")) {
      if (verbose) { message("Outer fold ", o, ": SRF ranking failed; skipping."); print(fit_rank) }
      next
    }
    
    if (is.null(fit_rank$importance) || length(fit_rank$importance) == 0) {
      if (verbose) message("Outer fold ", o, ": SRF returned empty importance; skipping.")
      next
    }
    
    vimp_sorted <- sort(fit_rank$importance, decreasing = TRUE)
    ranked_vars <- names(vimp_sorted)
    
    ranked_vars_outer[[o]] <- ranked_vars
    vimp_outer[[o]]        <- vimp_sorted
    
    if (verbose) {
      top_show <- min(10, length(ranked_vars))
      message("Outer fold ", o, ": top ranked vars: ",
              paste(ranked_vars[seq_len(top_show)], collapse = ", "))
    }
    
    # -----------------------------
    # (2) INNER CV (tune top_n)
    # -----------------------------
    if (verbose) message("Outer fold ", o, ": inner CV tuning top_n over {", paste(grid_top_n, collapse = ", "), "} …")
    
    folds_inner <- make_stratified_folds(
      dtrain_outer[[event_col]], K_inner, seed + o
    )
    
    grid_perf <- data.frame(
      top_n = grid_top_n,
      mean_bs = NA_real_,
      n_ok = 0L
    )
    
    for (g in seq_along(grid_top_n)) {
      top_n_g <- grid_top_n[g]
      top_vars_g <- ranked_vars[seq_len(min(top_n_g, length(ranked_vars)))]
      
      bs_vals <- rep(NA_real_, K_inner)
      ok      <- rep(FALSE,   K_inner)
      
      for (k in 1:K_inner) {
        tr_idx <- which(folds_inner != k)
        va_idx <- which(folds_inner == k)
        
        dtr <- dtrain_outer[tr_idx, , drop = FALSE]
        dva <- dtrain_outer[va_idx, , drop = FALSE]
        
        ev_dtr <- sum(dtr[[event_col]] == 1, na.rm = TRUE)
        ev_dva <- sum(dva[[event_col]] == 1, na.rm = TRUE)
        if (ev_dtr < 2 || ev_dva < 1) next
        
        max_va_time <- max(dva[[time_col]], na.rm = TRUE)
        times_inner <- times[times <= max_va_time]
        if (length(times_inner) == 0) next
        
        safe_top <- paste0("`", top_vars_g, "`")
        f_cox <- stats::as.formula(
          paste0("Surv(", time_col, ", ", event_col, ") ~ ",
                 paste(safe_top, collapse = " + "))
        )
        
        fit_cox <- try(coxph(f_cox, data = dtr, x = TRUE), silent = TRUE)
        if (inherits(fit_cox, "try-error")) next
        
        dscore <- dva
        dscore$time  <- dscore[[time_col]]
        dscore$event <- dscore[[event_col]]
        
        sc <- try(
          riskRegression::Score(
            object = list("HybridCox" = fit_cox),
            formula = surv_formula_score,
            data = dscore,
            metrics = "Brier",
            times = times_inner,
            null.model = TRUE,
            summary = "IPA"
          ),
          silent = TRUE
        )
        if (inherits(sc, "try-error")) next
        
        bs_model <- sc$Brier$score$Brier[sc$Brier$score$model == "HybridCox"]
        if (length(bs_model) == 0 || all(!is.finite(bs_model))) next
        
        bs_vals[k] <- mean(bs_model, na.rm = TRUE)
        ok[k] <- TRUE
      }
      
      grid_perf$mean_bs[g] <- mean(bs_vals, na.rm = TRUE)
      grid_perf$n_ok[g]    <- sum(ok)
      
      if (verbose) {
        message("  inner: top_n=", top_n_g,
                " | mean_BS=", round(grid_perf$mean_bs[g], 4),
                " | valid_folds=", grid_perf$n_ok[g], "/", K_inner)
      }
    }
    
    inner_grid_tables[[o]] <- grid_perf
    
    ok_rows <- grid_perf[grid_perf$n_ok > 0 & is.finite(grid_perf$mean_bs), ]
    if (nrow(ok_rows) == 0) {
      if (verbose) message("Outer fold ", o, ": inner tuning failed (no valid folds); skipping.")
      next
    }
    
    best_row <- ok_rows[which.min(ok_rows$mean_bs), ]
    best_topn_outer[o] <- best_row$top_n
    
    top_vars_best <- ranked_vars[
      seq_len(min(best_row$top_n, length(ranked_vars)))
    ]
    top_vars_outer[[o]] <- top_vars_best
    
    if (verbose) {
      message("Outer fold ", o, ": best_top_n=", best_row$top_n,
              " (inner mean_BS=", round(best_row$mean_bs, 4),
              ", valid_folds=", best_row$n_ok, "/", K_inner, ")")
    }
    
    # -----------------------------
    # (3) Refit Cox on full outer-train
    # -----------------------------
    safe_best <- paste0("`", top_vars_best, "`")
    f_best <- stats::as.formula(
      paste0("Surv(", time_col, ", ", event_col, ") ~ ",
             paste(safe_best, collapse = " + "))
    )
    
    fit_outer <- try(coxph(f_best, data = dtrain_outer, x = TRUE), silent = TRUE)
    if (inherits(fit_outer, "try-error")) {
      if (verbose) { message("Outer fold ", o, ": Cox refit failed; skipping."); print(fit_outer) }
      next
    }
    
    # -----------------------------
    # (4) Outer test evaluation
    # -----------------------------
    dscore_test <- dtest_outer
    dscore_test$time  <- dscore_test[[time_col]]
    dscore_test$event <- dscore_test[[event_col]]
    
    sc_out <- try(
      riskRegression::Score(
        object = list("HybridCox" = fit_outer),
        formula = surv_formula_score,
        data = dscore_test,
        metrics = "Brier",
        times = times_outer,
        null.model = TRUE,
        summary = "IPA"
      ),
      silent = TRUE
    )
    if (inherits(sc_out, "try-error")) {
      if (verbose) { message("Outer fold ", o, ": Score failed; skipping."); print(sc_out) }
      next
    }
    
    bs_model <- sc_out$Brier$score$Brier[sc_out$Brier$score$model == "HybridCox"]
    bs_null  <- sc_out$Brier$score$Brier[sc_out$Brier$score$model == "Null model"]
    
    bs_outer[o]      <- mean(bs_model, na.rm = TRUE)
    bs_null_outer[o] <- mean(bs_null,  na.rm = TRUE)
    ipa_outer[o]     <- 1 - (bs_outer[o] / bs_null_outer[o])
    
    if (verbose) {
      message("OUTER fold ", o,
              " | best_top_n=", best_topn_outer[o],
              " | BS=", round(bs_outer[o], 4),
              " | BS_null=", round(bs_null_outer[o], 4),
              " | IPA=", round(ipa_outer[o], 4))
    }
  }
  
  if (verbose) {
    message("\nScored folds: ", sum(is.finite(bs_outer)), "/", K_outer)
  }
  
  list(
    summary = data.frame(
      BS_mean = mean(bs_outer, na.rm = TRUE),
      BS_sd   = stats::sd(bs_outer, na.rm = TRUE),
      BS_null_mean = mean(bs_null_outer, na.rm= TRUE),
      BS_null_sd   = stats::sd(bs_null_outer, na.rm= TRUE),
      IPA_mean = mean(ipa_outer, na.rm = TRUE),
      IPA_sd   = stats::sd(ipa_outer, na.rm = TRUE),
      folds_scored = sum(is.finite(bs_outer)),
      stringsAsFactors = FALSE
    ),
    fold_results = data.frame(
      fold = 1:K_outer,
      BS = bs_outer,
      BS_null = bs_null_outer,
      IPA = ipa_outer,
      best_top_n = best_topn_outer,
      stringsAsFactors = FALSE
    ),
    ranked_vars_outer = ranked_vars_outer,
    vimp_outer = vimp_outer,
    top_vars_outer = top_vars_outer,
    inner_grid_tables = inner_grid_tables
  )
}



res_hybrid_clin <- doublecv_hybrid_srf_topn_cox_brier(
  df = df_imputed,
  vars = clinical_vars,
  time_col = "tfs_months",
  event_col = "event_death",
  K_outer = 10,
  K_inner = 5,
  times = seq(12,60,12),
  grid_top_n = c(5,10,15,20),
  importance_rank = "anti",   # fast
  ntree_rank = 500,
  seed = 123,
  verbose = TRUE
)

res_hybrid_clin$summary
head(res_hybrid_clin$fold_results, 10)

res_hybrid_clin$top_vars_outer[[1]]


# ============================================================
# 1) Runner: Hybrid SRF→Cox double-CV over var-group combos
#    (calls the engine doublecv_hybrid_srf_topn_cox_brier)
# ============================================================
run_doublecv_hybrid_over_combos <- function(
    df,
    var_groups,
    time_col    = "tfs_months",
    event_col   = "event_death",
    combo_sizes = 1:length(var_groups),
    
    # pass-through to engine
    K_outer = 10,
    K_inner = 5,                 # keep fixed at 5
    times   = seq(12, 60, 12),
    grid_top_n = c(5, 10, 15, 20),
    
    # SRF ranking (outer only)
    ntree_rank = 500,
    mtry       = NULL,
    nodesize   = 15,
    nsplit     = 10,
    importance_rank = "anti",
    
    # filters (like other runners)
    min_rows   = 30,
    min_events = 10,
    
    seed    = 123,
    verbose = TRUE
) {
  combos <- make_var_group_combos(var_groups, sizes = combo_sizes)
  
  summary_rows <- list()
  details_rows <- list()
  
  combo_id <- 1L
  
  for (groups_in in combos) {
    
    combo_label <- paste(groups_in, collapse = "+")
    vars <- unique(unlist(var_groups[groups_in], use.names = FALSE))
    
    if (verbose) {
      message("\n==============================")
      message("COMBO ", combo_id, ": ", combo_label, " (n_vars=", length(vars), ")")
      message("==============================")
    }
    
    # check columns exist
    needed <- c(time_col, event_col, vars)
    missing_cols <- setdiff(needed, colnames(df))
    if (length(missing_cols) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        combo_id = combo_id,
        groups   = combo_label,
        n_vars   = length(vars),
        rows_used = NA_integer_,
        events    = NA_integer_,
        folds_scored = 0L,
        BS_mean = NA_real_, BS_sd = NA_real_,
        BS_null_mean = NA_real_, BS_null_sd = NA_real_,
        IPA_mean = NA_real_, IPA_sd = NA_real_,
        status = paste("missing:", paste(missing_cols, collapse = ",")),
        stringsAsFactors = FALSE
      )
      details_rows[[as.character(combo_id)]] <- list(
        groups = groups_in,
        vars   = vars,
        result = NULL
      )
      combo_id <- combo_id + 1L
      next
    }
    
    # combo-specific complete cases (for min_rows/min_events check only)
    dsub_check <- df[stats::complete.cases(df[, needed, drop = FALSE]), needed, drop = FALSE]
    n_rows   <- nrow(dsub_check)
    n_events <- sum(dsub_check[[event_col]] == 1, na.rm = TRUE)
    
    if (n_rows < min_rows || n_events < min_events) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        combo_id = combo_id,
        groups   = combo_label,
        n_vars   = length(vars),
        rows_used = n_rows,
        events    = n_events,
        folds_scored = 0L,
        BS_mean = NA_real_, BS_sd = NA_real_,
        BS_null_mean = NA_real_, BS_null_sd = NA_real_,
        IPA_mean = NA_real_, IPA_sd = NA_real_,
        status = "insufficient rows/events",
        stringsAsFactors = FALSE
      )
      details_rows[[as.character(combo_id)]] <- list(
        groups = groups_in,
        vars   = vars,
        result = NULL
      )
      combo_id <- combo_id + 1L
      next
    }
    
    # ---- RUN THE ENGINE ----
    res <- doublecv_hybrid_srf_topn_cox_brier(
      df        = df,
      vars      = vars,
      time_col  = time_col,
      event_col = event_col,
      K_outer   = K_outer,
      K_inner   = K_inner,
      times     = times,
      grid_top_n = grid_top_n,
      ntree_rank = ntree_rank,
      mtry       = mtry,
      nodesize   = nodesize,
      nsplit     = nsplit,
      importance_rank = importance_rank,
      seed      = seed,
      verbose   = verbose
    )
    
    # ---- summarize like RSF runner ----
    fr <- res$fold_results
    folds_scored <- sum(is.finite(fr$BS))
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      combo_id = combo_id,
      groups   = combo_label,
      n_vars   = length(vars),
      rows_used = n_rows,
      events    = n_events,
      folds_scored = folds_scored,
      BS_mean = mean(fr$BS, na.rm = TRUE),
      BS_sd   = stats::sd(fr$BS, na.rm = TRUE),
      BS_null_mean = mean(fr$BS_null, na.rm = TRUE),
      BS_null_sd   = stats::sd(fr$BS_null, na.rm = TRUE),
      IPA_mean = mean(fr$IPA, na.rm = TRUE),
      IPA_sd   = stats::sd(fr$IPA, na.rm = TRUE),
      status   = "ok",
      stringsAsFactors = FALSE
    )
    
    details_rows[[as.character(combo_id)]] <- list(
      groups = groups_in,
      vars   = vars,
      result = res
    )
    
    combo_id <- combo_id + 1L
  }
  
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  
  list(
    summary = summary_df,
    details = details_rows
  )
}

var_groups <- list(
  lab      = lab_vars,
  ba       = ba_vars,
  ram      = ram_vars,
  clinical = clinical_vars
)

res_hybrid_all <- run_doublecv_hybrid_over_combos(
  df = df_imputed,
  var_groups = var_groups,
  combo_sizes = 1:4,
  K_outer = 10,
  K_inner = 5,
  times = seq(12, 60, 12),
  grid_top_n = c(5, 10, 15, 20),
  min_rows = 150,
  min_events = 10,
  importance_rank = "anti",
  seed = 123,
  verbose = TRUE
)

res_hybrid_all$summary


## -------------------------
## MULTI-PANEL SCATTERPLOT SETUP (2x3)
## -------------------------

par(
  mfrow = c(2, 3),   # draw 6 plots on one page (2 rows x 3 columns)
  pty = "s",         # force square plotting region (helps comparability across panels)
  mar = c(5, 5, 2, 1),   # margins (bottom, left, top, right)
  cex.axis = 1.2,        # axis tick label size
  cex.lab  = 1.4,        # axis title size
  font.axis=2,           # bold tick labels
  font.lab=2             # bold axis titles
)

## -------------------------
## PAIRWISE METHOD COMPARISONS (IPA mean) — SCATTERPLOTS
## Each point = one predictor-group combination (combo_id).
## Diagonal line y=x is "equal performance"; above line => y-method better; below => x-method better.
## asp=1 keeps x/y scales identical.
## -------------------------

plot(
  res_srf_all[["summary"]][["IPA_mean"]],
  res_ens2_all[["summary"]][["IPA_mean"]],
  xlim = c(-0.15,0.15), ylim = c(-0.15,0.15),  # fixed axis range for fair comparison
  xlab = "Mean IPA (RSF)",
  ylab = "Mean IPA (Tree-Cluster Cox)",
  pch = 19, cex = 1.6,  # point style and size
  asp = 1,              # 1:1 aspect ratio (so y=x is 45 degrees visually)
  mtext("(A)", side = 3, line = 0.3, adj = 0, font = 2)  # panel label
)
abline(0,1)  # reference line y = x


plot(
  res_srf_all[["summary"]][["IPA_mean"]],
  results_doublecv_imputed[["summary"]][["IPA_mean"]],
  xlim = c(-0.15,0.15), ylim = c(-0.15,0.15),
  xlab = "Mean IPA (RSF)",
  ylab = "Mean IPA (LASSO-Cox)",
  pch = 19, cex = 1.6,
  asp = 1,
  mtext("(B)", side = 3, line = 0.3, adj = 0, font = 2)
)
abline(0,1)


plot(
  res_srf_all[["summary"]][["IPA_mean"]],
  res_hybrid_all[["summary"]][["IPA_mean"]],
  xlim = c(-0.15,0.15), ylim = c(-0.15,0.15),
  xlab = "Mean IPA (RSF)",
  ylab = "Mean IPA (Hybrid RSF-Cox )",
  pch = 19, cex = 1.6,
  asp = 1,
  mtext("(C)", side = 3, line = 0.3, adj = 0, font = 2)
)
abline(0,1)


plot(
  res_ens2_all[["summary"]][["IPA_mean"]],
  results_doublecv_imputed[["summary"]][["IPA_mean"]],
  xlim =c(-0.15,0.15), ylim = c(-0.15,0.15),
  xlab = "Mean IPA (Tree-Cluster Cox)",
  ylab = "Mean IPA (LASSO-Cox)",
  pch = 19, cex = 1.6,
  asp = 1,
  mtext("(D)", side = 3, line = 0.3, adj = 0, font = 2)
)
abline(0,1)


plot(
  res_ens2_all[["summary"]][["IPA_mean"]],
  res_hybrid_all[["summary"]][["IPA_mean"]],
  xlim = c(-0.15,0.15), ylim = c(-0.15,0.15),
  xlab = "Mean IPA (Tree-Cluster Cox )",
  ylab = "Mean IPA (Hybrid RSF-Cox)",
  pch = 19, cex = 1.6,
  asp = 1,
  mtext("(E)", side = 3, line = 0.3, adj = 0, font = 2)
)
abline(0,1)


plot(
  results_doublecv_imputed[["summary"]][["IPA_mean"]],
  res_hybrid_all[["summary"]][["IPA_mean"]],
  xlim = c(-0.15,0.15), ylim = c(-0.15,0.15),
  xlab = "Mean IPA (LASSO-Cox)",
  ylab = "Mean IPA (Hybrid RSF-Cox)",
  pch = 19, cex = 1.6,
  asp = 1,
  mtext("(F)", side = 3, line = 0.3, adj = 0, font = 2)
)
abline(0,1)


## -------------------------
## RESET PLOTTING SETTINGS BACK TO DEFAULTS
## -------------------------
par(
  mfrow = c(1, 1),  # back to single plot per page
  cex.axis = 1,
  cex.lab = 1,
  font.axis = 1,
  font.lab = 1
)


## -------------------------
## RESHAPE RESULTS INTO "LONG" FORMAT (tidy) FOR COMPARISON
## -------------------------

to_long <- function(res, method_name) {
  # Extracts a small standardized table from each result object:
  # combo_id + groups identify the predictor-set, IPA_mean is the metric, method tags the source.
  res$summary %>%
    select(combo_id, groups, IPA_mean) %>%
    mutate(method = method_name)
}

# Stack all methods into one long table: one row per (combo_id, method).
ipa_long <- bind_rows(
  to_long(results_doublecv_imputed,  "LASSO"),
  to_long(res_srf_all,    "RSF"),
  to_long(res_ens2_all,   "Ensemble2"),
  to_long(res_hybrid_all, "Hybrid")
)


## -------------------------
## WIDEN INTO ONE ROW PER COMBINATION (15 x 4 TABLE)
## -------------------------

ipa_15x4 <- ipa_long %>%
  pivot_wider(
    id_cols    = c(combo_id, groups),      # keys for each predictor-group combination
    names_from = method,                   # create one column per method
    values_from = IPA_mean                 # fill with IPA_mean values
  ) %>%
  arrange(combo_id)                        # sort by combination id

## -------------------------
## GLOBAL BEST IPA (ACROSS ALL METHODS AND COMBINATIONS)
## -------------------------

global_max_info <- ipa_15x4 %>%
  pivot_longer(
    cols = c(LASSO, RSF, Ensemble2, Hybrid),  # reshape back to method/value columns
    names_to = "method",
    values_to = "IPA"
  ) %>%
  filter(IPA == max(IPA, na.rm = TRUE))       # keep row(s) achieving global maximum

global_max_info  # shows which combo_id + method achieved the best overall IPA


## -------------------------
## ROW-WISE BEST METHOD PER COMBINATION
## -------------------------

rowwise_best <- ipa_15x4 %>%
  rowwise() %>%  # treat each row independently
  mutate(
    best_IPA = max(c(LASSO, RSF, Ensemble2, Hybrid), na.rm = TRUE),  # best value in that row
    best_method = c("LASSO","RSF","Ensemble2","Hybrid")[
      which.max(c(LASSO, RSF, Ensemble2, Hybrid))  # method that attains the max (ties pick first)
    ]
  ) %>%
  ungroup()

rowwise_best  # per combo_id: best method and its IPA


## -------------------------
## COLUMN-WISE BEST IPA PER METHOD (METHOD-SPECIFIC MAX OVER COMBINATIONS)
## -------------------------

colwise_best <- ipa_15x4 %>%
  summarise(
    LASSO_max      = max(LASSO, na.rm = TRUE),
    RSF_max        = max(RSF, na.rm = TRUE),
    Ensemble2_max  = max(Ensemble2, na.rm = TRUE),
    Hybrid_max     = max(Hybrid, na.rm = TRUE)
  )

colwise_best  # gives each method’s best achieved IPA across all predictor combinations


## -------------------------
## FOR EACH METHOD: WHICH COMBINATION ACHIEVED ITS MAX?
## -------------------------

colwise_best_info <- ipa_15x4 %>%
  pivot_longer(
    cols = c(LASSO, RSF, Ensemble2, Hybrid),
    names_to = "method",
    values_to = "IPA"
  ) %>%
  group_by(method) %>%
  filter(IPA == max(IPA, na.rm = TRUE))  # keep the top-performing combo(s) for each method
%>%
  ungroup()

colwise_best_info  # returns the combo_id(s) that maximize IPA for each method


## -------------------------
## INSPECT SUMMARY TABLES (DEBUG / REPORTING)
## -------------------------

results_doublecv_imputed$summary  # LASSO-Cox summary metrics per combination
res_srf_all$summary               # RSF summary metrics per combination
res_ens2_all$summary              # Ensemble2 summary metrics per combination
res_hybrid_all$summary            # Hybrid summary metrics per combination


## -------------------------
## INSPECT SELECTED VARIABLES AND LAMBDA DISTRIBUTIONS (LASSO DETAILS)
## -------------------------

# Shows which predictors were selected in each CV fold for combination "6" and "7"
# (useful to assess selection stability).
results_doublecv_imputed$details[["6"]]$selected_vars_per_fold
results_doublecv_imputed$details[["7"]]$selected_vars_per_fold

# Summarises the distribution of the chosen lambda.min across folds
# (helps see how variable/unstable the tuning is).
summary(results_doublecv_imputed$details[["6"]]$lambda_min)
summary(results_doublecv_imputed$details[["7"]]$lambda_min)

# vimp list from combo 7
vimp_list <- res_srf_all[["details"]][["7"]][["result"]][["vimp_outer"]]

K_outer <- length(vimp_list)   # should be 10
top_k   <- 10                  # top-10 per fold

# 1) Get top-10 variable names per fold (by VIMP, descending)
top_vars_per_fold <- lapply(vimp_list, function(v) {
  v <- v[is.finite(v)]
  names(sort(v, decreasing = TRUE))[seq_len(min(top_k, length(v)))]
})

# 2) Count how often each variable appears in top-10 across folds
tab <- sort(table(unlist(top_vars_per_fold)), decreasing = TRUE)

stability_df <- data.frame(
  variable = names(tab),
  n_folds_top10 = as.integer(tab),
  prop_top10 = as.integer(tab) / K_outer
)

# 3) Keep variables that appear in top-10 in >=70% of folds
stable_70 <- subset(stability_df, prop_top10 >= 0.70)

stable_70

stable_70$percent_top10 <- round(100 * stable_70$prop_top10, 0)
stable_70 <- stable_70[order(-stable_70$n_folds_top10, stable_70$variable), ]

stable_70



# vimp list from combo 13
vimp_list13<- res_srf_all[["details"]][["13"]][["result"]][["vimp_outer"]]

K_outer <- length(vimp_list13)   # should be 10
top_k   <- 10                  # top-10 per fold

# 1) Get top-10 variable names per fold (by VIMP, descending)
top_vars_per_fold <- lapply(vimp_list13, function(v) {
  v <- v[is.finite(v)]
  names(sort(v, decreasing = TRUE))[seq_len(min(top_k, length(v)))]
})

# 2) Count how often each variable appears in top-10 across folds
tab <- sort(table(unlist(top_vars_per_fold)), decreasing = TRUE)

stability_df <- data.frame(
  variable = names(tab),
  n_folds_top10 = as.integer(tab),
  prop_top10 = as.integer(tab) / K_outer
)

# 3) Keep variables that appear in top-10 in >=70% of folds
stable_70 <- subset(stability_df, prop_top10 >= 0.70)

stable_70

stable_70$percent_top10 <- round(100 * stable_70$prop_top10, 0)
stable_70 <- stable_70[order(-stable_70$n_folds_top10, stable_70$variable), ]

stable_70


# vimp list from combo 6
vimp_list6<- res_srf_all[["details"]][["6"]][["result"]][["vimp_outer"]]

K_outer <- length(vimp_list6)   # should be 10
top_k   <- 10                  # top-10 per fold

# 1) Get top-10 variable names per fold (by VIMP, descending)
top_vars_per_fold <- lapply(vimp_list6, function(v) {
  v <- v[is.finite(v)]
  names(sort(v, decreasing = TRUE))[seq_len(min(top_k, length(v)))]
})

# 2) Count how often each variable appears in top-10 across folds
tab <- sort(table(unlist(top_vars_per_fold)), decreasing = TRUE)

stability_df <- data.frame(
  variable = names(tab),
  n_folds_top10 = as.integer(tab),
  prop_top10 = as.integer(tab) / K_outer
)

# 3) Keep variables that appear in top-10 in >=70% of folds
stable_70 <- subset(stability_df, prop_top10 >= 0.70)

stable_70

stable_70$percent_top10 <- round(100 * stable_70$prop_top10, 0)
stable_70 <- stable_70[order(-stable_70$n_folds_top10, stable_70$variable), ]

stable_70


library(randomForestSRC)
library(rpart)
library(survival)

# ---- choose which combo to illustrate ----
combo_id <- 2  #clinical (example)
vars_illustrate <- res_ens2_all$details[[as.character(combo_id)]]$vars

time_col  <- "tfs_months"
event_col <- "event_death"

# ---- make complete-case subset on needed cols (same as the engine) ----
needed <- c(time_col, event_col, vars_illustrate)
dsub_full <- df_imputed[complete.cases(df_imputed[, needed, drop = FALSE]), needed, drop = FALSE]
dsub_full[[time_col]] <- ifelse(dsub_full[[time_col]] <= 0, 1e-6, dsub_full[[time_col]])

# ---- 1) RSF ranking on full data (permute importance) ----
set.seed(123)
mtry_rank <- max(1, floor(sqrt(length(vars_illustrate))))

f_rsf <- as.formula(
  paste0("Surv(", time_col, ", ", event_col, ") ~ ",
         paste0("`", vars_illustrate, "`", collapse = " + "))
)

fit_rank_full <- rfsrc(
  formula    = f_rsf,
  data       = dsub_full,
  ntree      = 500,
  mtry       = mtry_rank,
  nodesize   = 15,
  nsplit     = 10,
  importance = "permute"
)

vimp_full <- sort(fit_rank_full$importance, decreasing = TRUE)

# ---- 2) pick top_n and tree depth for the illustration ----
# Option A (simple): choose top_n=10, depth=1 (matches most of the outer folds)
top_n_ill <- 10
depth_ill <- 2

top_vars_ill <- names(vimp_full)[seq_len(min(top_n_ill, length(vimp_full)))]

# ---- 3) Stage 2 fit: tree clustering + Cox-per-cluster ----
st2_full <- build_ens2_stage2_tree_clustercox(
  df = dsub_full,
  top_vars = top_vars_ill,
  time_col = time_col,
  event_col = event_col,
  tree_maxdepth  = depth_ill,
  tree_minsplit  = 30,
  tree_minbucket = 15,
  tree_cp        = 0.01,
  min_cluster_events = 5,
  seed = 123
)

# ---- 4) show what the clusters look like ----
st2_full$tree$cluster_table
st2_full$cox_by_cluster$dropped_clusters
st2_full$top_vars

pdf("ens2_tree_example.pdf", width = 6, height = 4)

par(mar = c(1, 1, 1, 1))

rpart.plot(
  st2_full$tree$fit,
  roundint = FALSE,
  extra = 102,
  under = TRUE,
  faclen = 0,
  compress = TRUE,
  branch = 0.4,
  tweak = 0.9
)

dev.off()
st2_full$cox_by_cluster$formula

# list available clusters with fitted Cox models
names(st2_full$cox_by_cluster$models)

# print model summaries cluster-by-cluster
for (cl in names(st2_full$cox_by_cluster$models)) {
  cat("\n========================\n")
  cat("Cluster:", cl, "\n")
  cat("========================\n")
  print(summary(st2_full$cox_by_cluster$models[[cl]]))
}


