library(T2DMTEST)

renv::restore()

# Maximum number of cores to be used:
maxCores <- parallel::detectCores()

# The folder where the study intermediate and result files will be written:
outputFolder <- ""

# Details for connecting to the server:
connectionDetails <- DatabaseConnector::createConnectionDetails(dbms = "",
                                                                server = "",
                                                                user = "",
                                                                password = "",
                                                                port = "",
                                                                pathToDriver = "")


# The name of the database schema where the CDM data can be found:
cdmDatabaseSchema <- ""
vocabulary_database_schema <- ""
# The name of the database schema and table where the study-specific cohorts will be instantiated:
cohortDatabaseSchema <- ""
cohortTable <- ""

# Some meta-information that will be used by the export function:
databaseId <- ""
databaseName <- ""
databaseDescription <- ""
# For Oracle: define a schema that can be used to emulate temp tables:
oracleTempSchema <- NULL

execute(connectionDetails = connectionDetails,
        cdmDatabaseSchema = cdmDatabaseSchema,
        cohortDatabaseSchema = cohortDatabaseSchema,
        cohortTable = cohortTable,
        oracleTempSchema = oracleTempSchema,
        outputFolder = outputFolder,
        databaseId = databaseId,
        databaseName = databaseName,
        databaseDescription = databaseDescription,
        createCohorts = TRUE,
        synthesizePositiveControls = TRUE,
        runAnalyses = TRUE,
        packageResults = TRUE,
        maxCores = maxCores)

#cohort method 데이터 추출(위치와 파일명 확인)
RDSpath <- paste(outputFolder,sep="","/cmOutput/Ps_l1_s1_p1_t293_c368_o253.rds")
df_cm <- readRDS(RDSpath)

#합병증 개수별 코호트 생성
conn <- DatabaseConnector::connect(connectionDetails)

source("extras/target.R")
source("extras/single.R")
source("extras/double.R")
source("extras/triple.R")

#합병증 개수별 코호트 환자 데이터 추출(추후 Characterization에 활용)
df_0 <- renderTranslateQuerySql(conn,
   'SELECT [cohort_definition_id]
      ,[subject_id],
	  person_id,year_of_birth,gender_concept_id
  FROM @target_database_schema.@target_cohort_table as a
    inner join (select *
	from @cdm_database_schema.PERSON) as b
	on a.subject_id = b.person_id
  where cohort_definition_id = 0',
  cdm_database_schema = cdmDatabaseSchema,
  target_cohort_table = cohortTable,
  target_database_schema = cohortDatabaseSchema)

df_1 <- renderTranslateQuerySql(conn,
                                'SELECT [cohort_definition_id]
      ,[subject_id],
	  person_id,year_of_birth,gender_concept_id
  FROM @target_database_schema.@target_cohort_table as a
    inner join (select *
	from @cdm_database_schema.PERSON) as b
	on a.subject_id = b.person_id
  where cohort_definition_id = 1',
  cdm_database_schema = cdmDatabaseSchema,
  target_cohort_table = cohortTable,
  target_database_schema = cohortDatabaseSchema)

df_2 <- renderTranslateQuerySql(conn,
                                'SELECT [cohort_definition_id]
      ,[subject_id],
	  person_id,year_of_birth,gender_concept_id
  FROM @target_database_schema.@target_cohort_table as a
    inner join (select *
	from @cdm_database_schema.PERSON) as b
	on a.subject_id = b.person_id
  where cohort_definition_id = 2',
  cdm_database_schema = cdmDatabaseSchema,
  target_cohort_table = cohortTable,
  target_database_schema = cohortDatabaseSchema)

df_3 <- renderTranslateQuerySql(conn,
                                'SELECT [cohort_definition_id]
      ,[subject_id],
	  person_id,year_of_birth,gender_concept_id
  FROM @target_database_schema.@target_cohort_table as a
    inner join (select *
	from @cdm_database_schema.PERSON) as b
	on a.subject_id = b.person_id
  where cohort_definition_id = 3',
  cdm_database_schema = cdmDatabaseSchema,
  target_cohort_table = cohortTable,
  target_database_schema = cohortDatabaseSchema)



#cohort method 환자군 분리(R의경우 반복문의 속도가 느린편이라 list 형식으로 변경후 비교)
df_treatment <- df_cm$treatment
df_id <- df_cm$subjectId
df_id_0 <- df_0$SUBJECT_ID
df_id_1 <- df_1$SUBJECT_ID
df_id_2 <- df_2$SUBJECT_ID
df_id_3 <- df_3$SUBJECT_ID

#합병증 개수 0개 환자
for(i in 1:length(df_id)){
  for( j in 1:length(df_id_0)){
    if(df_id[i] == df_id_0[j])
      df_treatment[i] <- 0
  }
}

#합병증 개수 1개 환자
for(i in 1:length(df_id)){
  for( j in 1:length(df_id_1)){
    if(df_id[i] == df_id_1[j])
      df_treatment[i] <- 1
  }
}

#합병증 개수 2개 환자
for(i in 1:length(df_id)){
  for( j in 1:length(df_id_2)){
    if(df_id[i] == df_id_2[j])
      df_treatment[i] <- 2
  }
}

#합병증 개수 3개 환자
for(i in 1:length(df_id)){
  for( j in 1:length(df_id_3)){
    if(df_id[i] == df_id_3[j])
      df_treatment[i] <- 3
  }
}

#합병증 개수 변경후 저장
df_cm$treatment <- df_treatment

write.csv(df_cm,file="df_cm.csv")


#합병증 개수별 코호트 Target, Comparator 분리
ple_0 <- df_cm %>% filter(df_cm$treatment == 0)
ple_0$treatment <- 1
ple_1 <- df_cm %>% filter(df_cm$treatment == 1)
ple_1$treatment <- 0
ple_2 <- df_cm %>% filter(df_cm$treatment == 2)
ple_2$treatment <- 0
ple_3 <- df_cm %>% filter(df_cm$treatment == 3)
ple_3$treatment <- 0

#합병증 개수별 코호트 조합
df_solo <- rbind(ple_0,ple_1)
df_double <- rbind(ple_0,ple_2)
df_triple <- rbind(ple_0,ple_3)


#코호트별 성향점수매칭
matchedPop_1 <- matchOnPs(population = df_solo, caliper = 0,
                          caliperScale = "standardized logit", maxRatio = 0)

matchedPop_2 <- matchOnPs(population = df_double, caliper = 0,
                          caliperScale = "standardized logit", maxRatio = 0)

matchedPop_3 <- matchOnPs(population = df_triple, caliper = 0,
                          caliperScale = "standardized logit", maxRatio = 0)


#COX 모델 생성
outcomeModel_1 <- fitOutcomeModel(population = df_solo,
                                  modelType = "cox",
                                  stratified = FALSE)
outcomeModel_2 <- fitOutcomeModel(population = df_double,
                                  modelType = "cox",
                                  stratified = FALSE)
outcomeModel_3 <- fitOutcomeModel(population = df_triple,
                                  modelType = "cox",
                                  stratified = FALSE)
#COX 모델 출력
outcomeModel_1
outcomeModel_2
outcomeModel_3

#AttritionDiagram(Original cohorts내용은 분리전 초기값 출력)
drawAttritionDiagram(matchedPop_1)
drawAttritionDiagram(matchedPop_2)
drawAttritionDiagram(matchedPop_3)


#prepareKaplanMeier
prepareKaplanMeier <- function(population) {
  dataCutoff <- 0.9
  population$y <- 0
  population$y[population$outcomeCount != 0] <- 1
  if (is.null(population$stratumId) || length(unique(population$stratumId)) == nrow(population)/2) {
    sv <- survival::survfit(survival::Surv(survivalTime, y) ~ treatment, population, conf.int = TRUE)
    idx <- summary(sv, censored = T)$strata == "treatment=1"
    survTarget <- tibble::tibble(time = sv$time[idx],
                                 targetSurvival = sv$surv[idx],
                                 targetSurvivalLb = sv$lower[idx],
                                 targetSurvivalUb = sv$upper[idx])
    idx <- summary(sv, censored = T)$strata == "treatment=0"
    survComparator <- tibble::tibble(time = sv$time[idx],
                                     comparatorSurvival = sv$surv[idx],
                                     comparatorSurvivalLb = sv$lower[idx],
                                     comparatorSurvivalUb = sv$upper[idx])
    data <- merge(survTarget, survComparator, all = TRUE)
  } else {
    population$stratumSizeT <- 1
    strataSizesT <- aggregate(stratumSizeT ~ stratumId, population[population$treatment == 1, ], sum)
    if (max(strataSizesT$stratumSizeT) == 1) {
      # variable ratio matching: use propensity score to compute IPTW
      if (is.null(population$propensityScore)) {
        stop("Variable ratio matching detected, but no propensity score found")
      }
      weights <- aggregate(propensityScore ~ stratumId, population, mean)
      if (max(weights$propensityScore) > 0.99999) {
        return(NULL)
      }
      weights$weight <- weights$propensityScore / (1 - weights$propensityScore)
    } else {
      # stratification: infer probability of treatment from subject counts
      strataSizesC <- aggregate(stratumSizeT ~ stratumId, population[population$treatment == 0, ], sum)
      colnames(strataSizesC)[2] <- "stratumSizeC"
      weights <- merge(strataSizesT, strataSizesC)
      if (nrow(weights) == 0) {
        warning("No shared strata between target and comparator")
        return(NULL)
      }
      weights$weight <- weights$stratumSizeT/weights$stratumSizeC
    }
    population <- merge(population, weights[, c("stratumId", "weight")])
    population$weight[population$treatment == 1] <- 1
    idx <- population$treatment == 1
    survTarget <- CohortMethod:::adjustedKm(weight = population$weight[idx],
                                            time = population$survivalTime[idx],
                                            y = population$y[idx])
    survTarget$targetSurvivalUb <- survTarget$s^exp(qnorm(0.975)/log(survTarget$s) * sqrt(survTarget$var)/survTarget$s)
    survTarget$targetSurvivalLb <- survTarget$s^exp(qnorm(0.025)/log(survTarget$s) * sqrt(survTarget$var)/survTarget$s)
    survTarget$targetSurvivalLb[survTarget$s > 0.9999] <- survTarget$s[survTarget$s > 0.9999]
    survTarget$targetSurvival <- survTarget$s
    survTarget$s <- NULL
    survTarget$var <- NULL
    idx <- population$treatment == 0
    survComparator <- CohortMethod:::adjustedKm(weight = population$weight[idx],
                                                time = population$survivalTime[idx],
                                                y = population$y[idx])
    survComparator$comparatorSurvivalUb <- survComparator$s^exp(qnorm(0.975)/log(survComparator$s) *
                                                                  sqrt(survComparator$var)/survComparator$s)
    survComparator$comparatorSurvivalLb <- survComparator$s^exp(qnorm(0.025)/log(survComparator$s) *
                                                                  sqrt(survComparator$var)/survComparator$s)
    survComparator$comparatorSurvivalLb[survComparator$s > 0.9999] <- survComparator$s[survComparator$s >
                                                                                         0.9999]
    survComparator$comparatorSurvival <- survComparator$s
    survComparator$s <- NULL
    survComparator$var <- NULL
    data <- merge(survTarget, survComparator, all = TRUE)
  }
  data <- data[, c("time", "targetSurvival", "targetSurvivalLb", "targetSurvivalUb", "comparatorSurvival", "comparatorSurvivalLb", "comparatorSurvivalUb")]
  cutoff <- quantile(population$survivalTime, dataCutoff)
  data <- data[data$time <= cutoff, ]
  if (cutoff <= 300) {
    xBreaks <- seq(0, cutoff, by = 50)
  } else if (cutoff <= 600) {
    xBreaks <- seq(0, cutoff, by = 100)
  } else {
    xBreaks <- seq(0, cutoff, by = 250)
  }
  
  targetAtRisk <- c()
  comparatorAtRisk <- c()
  for (xBreak in xBreaks) {
    targetAtRisk <- c(targetAtRisk,
                      sum(population$treatment == 1 & population$survivalTime >= xBreak))
    comparatorAtRisk <- c(comparatorAtRisk,
                          sum(population$treatment == 0 & population$survivalTime >=
                                xBreak))
  }
  data <- merge(data, tibble::tibble(time = xBreaks,
                                     targetAtRisk = targetAtRisk,
                                     comparatorAtRisk = comparatorAtRisk), all = TRUE)
  if (is.na(data$targetSurvival[1])) {
    data$targetSurvival[1] <- 1
    data$targetSurvivalUb[1] <- 1
    data$targetSurvivalLb[1] <- 1
  }
  if (is.na(data$comparatorSurvival[1])) {
    data$comparatorSurvival[1] <- 1
    data$comparatorSurvivalUb[1] <- 1
    data$comparatorSurvivalLb[1] <- 1
  }
  idx <- which(is.na(data$targetSurvival))
  while (length(idx) > 0) {
    data$targetSurvival[idx] <- data$targetSurvival[idx - 1]
    data$targetSurvivalLb[idx] <- data$targetSurvivalLb[idx - 1]
    data$targetSurvivalUb[idx] <- data$targetSurvivalUb[idx - 1]
    idx <- which(is.na(data$targetSurvival))
  }
  idx <- which(is.na(data$comparatorSurvival))
  while (length(idx) > 0) {
    data$comparatorSurvival[idx] <- data$comparatorSurvival[idx - 1]
    data$comparatorSurvivalLb[idx] <- data$comparatorSurvivalLb[idx - 1]
    data$comparatorSurvivalUb[idx] <- data$comparatorSurvivalUb[idx - 1]
    idx <- which(is.na(data$comparatorSurvival))
  }
  data$targetSurvival <- round(data$targetSurvival, 4)
  data$targetSurvivalLb <- round(data$targetSurvivalLb, 4)
  data$targetSurvivalUb <- round(data$targetSurvivalUb, 4)
  data$comparatorSurvival <- round(data$comparatorSurvival, 4)
  data$comparatorSurvivalLb <- round(data$comparatorSurvivalLb, 4)
  data$comparatorSurvivalUb <- round(data$comparatorSurvivalUb, 4)
  
  # Remove duplicate (except time) entries:
  data <- data[order(data$time), ]
  data <- data[!duplicated(data[, -1]), ]
  return(data)
}

#plotKaplanMeier(y축 조정시 ylims 값 수정)
plotKaplanMeier <- function(kaplanMeier, targetName, comparatorName) {
  data <- rbind(data.frame(time = kaplanMeier$time,
                           s = kaplanMeier$targetSurvival,
                           lower = kaplanMeier$targetSurvivalLb,
                           upper = kaplanMeier$targetSurvivalUb,
                           strata = paste0(" ", targetName, "    ")),
                data.frame(time = kaplanMeier$time,
                           s = kaplanMeier$comparatorSurvival,
                           lower = kaplanMeier$comparatorSurvivalLb,
                           upper = kaplanMeier$comparatorSurvivalUb,
                           strata = paste0(" ", comparatorName)))
  
  xlims <- c(-max(data$time)/40, max(data$time))
  ylims <- c(0.7, 1)
  xLabel <- "Time in days"
  yLabel <- "Survival probability"
  xBreaks <- kaplanMeier$time[!is.na(kaplanMeier$targetAtRisk)]
  plot <- ggplot2::ggplot(data, ggplot2::aes(x = time,
                                             y = s,
                                             color = strata,
                                             fill = strata,
                                             ymin = lower,
                                             ymax = upper)) +
    ggplot2::geom_ribbon(color = rgb(0, 0, 0, alpha = 0)) +
    ggplot2::geom_step(size = 1) +
    ggplot2::scale_color_manual(values = c(rgb(0.8, 0, 0, alpha = 0.8),
                                           rgb(0, 0, 0.8, alpha = 0.8))) +
    ggplot2::scale_fill_manual(values = c(rgb(0.8, 0, 0, alpha = 0.3),
                                          rgb(0, 0, 0.8, alpha = 0.3))) +
    ggplot2::scale_x_continuous(xLabel, limits = xlims, breaks = xBreaks) +
    ggplot2::scale_y_continuous(yLabel, limits = ylims) +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
                   legend.position = "top",
                   legend.key.size = ggplot2::unit(1, "lines"),
                   plot.title = ggplot2::element_text(hjust = 0.5)) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(vjust = -10))
  
  targetAtRisk <- kaplanMeier$targetAtRisk[!is.na(kaplanMeier$targetAtRisk)]
  comparatorAtRisk <- kaplanMeier$comparatorAtRisk[!is.na(kaplanMeier$comparatorAtRisk)]
  labels <- data.frame(x = c(0, xBreaks, xBreaks),
                       y = as.factor(c("Number at risk",
                                       rep(targetName, length(xBreaks)),
                                       rep(comparatorName, length(xBreaks)))),
                       label = c("",
                                 formatC(targetAtRisk, big.mark = ",", mode = "integer"),
                                 formatC(comparatorAtRisk, big.mark = ",", mode = "integer")))
  labels$y <- factor(labels$y, levels = c(comparatorName, targetName, "Number at risk"))
  dataTable <- ggplot2::ggplot(labels, ggplot2::aes(x = x, y = y, label = label)) + ggplot2::geom_text(size = 3.5, vjust = 0.5) + ggplot2::scale_x_continuous(xLabel,
                                                                                                                                                              limits = xlims,
                                                                                                                                                              breaks = xBreaks) + ggplot2::theme(panel.grid.major = ggplot2::element_blank(),                                                                                                                                                                                  axis.title.y = ggplot2::element_blank(),
                                                                                                                                                                                                 axis.ticks = ggplot2::element_line(color = "white"))
  plots <- list(plot, dataTable)
  grobs <- widths <- list()
  for (i in 1:length(plots)) {
    grobs[[i]] <- ggplot2::ggplotGrob(plots[[i]])
    widths[[i]] <- grobs[[i]]$widths[2:5]
  }
  maxwidth <- do.call(grid::unit.pmax, widths)
  for (i in 1:length(grobs)) {
    grobs[[i]]$widths[2:5] <- as.list(maxwidth)
  }
  plot <- gridExtra::grid.arrange(grobs[[1]], grobs[[2]], heights = c(400, 100))
  return(plot)
}

#KaplanMeier 생존분석
pk1 <- prepareKaplanMeier(matchedPop_1)
pk2 <- prepareKaplanMeier(matchedPop_2)
pk3 <- prepareKaplanMeier(matchedPop_3)

plotKaplanMeier(pk1, "target","comparator")
plotKaplanMeier(pk2, "target","comparator")
plotKaplanMeier(pk3, "target","comparator")

#DB연결 종료
dbDisconnect(conn)

