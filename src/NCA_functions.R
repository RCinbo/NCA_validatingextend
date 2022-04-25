library(assertthat)
library(caret) # Confusion matrix maken
cleanmapdata <- function(data = data, points_id, tbltrans){
  map <- terra::extract(x = data,
                 y = terra::vect(
                   st_as_sf(x = points_id[, c('POINT_X','POINT_Y')],
                            coords = c("POINT_X", "POINT_Y"),
                            crs = "EPSG:31370"))
  )  #extract de waardes van de referentiedata
  names(map)[2] <- "landgebruik"
  map <- map %>%
    mutate(x = points_id[, 'POINT_X'],
           y = points_id[, 'POINT_X'],
           code = as.factor(landgebruik),
           landgebruik = recode_factor(code, "1" = "Open natuur", "2" = "Bos",
                                       "3" = "Grasland", "4" = "Akker",  "5" =
                                         "Urbaan", "6" = "Laag groen", "7" =
                                         "Hoog groen", "8" = "Water", "9" = "Overig"
           )) %>%
    left_join(tbltrans[,-1], by = c("code" = "lucode")) %>% droplevels()
  return(map)
}
#mapdata <- nara13$valid
#refdata <- as.factor(points_id$lu13oord)
#both arrays need to be factor variables with the same levels.
CalculateAccuracy <- function(mapdata, refdata){
  assert_that(length(mapdata) == length(refdata),
              msg = "Length of the arrays is not equal")
  assert_that(nlevels(mapdata) == nlevels(refdata) &
                all(levels(mapdata) %in% levels(refdata)),
              msg = "The data are not factors or don't have the same levels.")
  conf <- confusionMatrix(data = mapdata,
                          reference = refdata)
  #overall accuracy
  #sum(diag(conf$table))/sum(conf$table)
  return(conf)
}

library(assertthat)
library(caret) # Confusion matrix maken
#mapdata <- nara13$valid
#refdata <- as.factor(points_id$lu13oord)
#both arrays need to be factor variables with the same levels.
CalculateAccuracyChange <- function(mapdata1, refdata1, mapdata2, refdata2){
  assert_that(length(mapdata1) == length(refdata1) &
                length(mapdata1) == length(refdata2) &
                length(mapdata2) == length(refdata2),
              msg = "Length of the arrays is not equal")
  assert_that(nlevels(mapdata1) == nlevels(refdata1) &
                nlevels(mapdata1) == nlevels(refdata2) &
                nlevels(mapdata2) == nlevels(refdata2) &
                all(levels(mapdata1) %in% levels(refdata1)) &
                all(levels(mapdata1) %in% levels(refdata2)) &
                all(levels(mapdata2) %in% levels(refdata2)),
              msg = "The data are not factors or do not have the same levels.")


}

validationData <- function(){
  punten <- read_vc("validatiepunten", root = sprintf("%s/data/",here())) # Gevalideerde punten
  combine <- read_vc("combine", root = sprintf("%s/data/",here())) # Combine van de validatieklassenkaart van 2013
  # en 2016 -> geeft de oppervlakte van de landgebruiksveranderingen en
  # van de stabiele klassen
  # Omzettingstabel landgebruiken
  lu <- c(
    "Open natuur", "Bos", "Grasland", "Akker", "Urbaan", "Laag groen",
    "Hoog groen", "Water", "Overig"
  )
  lucode <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
  valid <- c(
    "Open natuur", "Hoog groen", "Open natuur", "Akker", "Urbaan",
    "Open natuur", "Hoog groen", "Water", "Overig"
  )
  valid_eng <- c(
    "Open nature", "High green", "Open nature", "Field", "Urban",
    "Open nature", "High green", "Water", "Other"
  )
  validcode <- c(1, 2, 1, 4, 5, 1, 2, 8, 9)
  tbltrans <- data.frame(lu, lucode, valid, valid_eng, validcode) %>%
    mutate(valid = as.factor(valid),
           valid_eng = as.factor(valid_eng),
           lucode = as.factor(lucode),
           validcode = as.factor(validcode),
           lu = as.factor(lu))

  # OPMAAK DATASET -----------------------------------------------------
  # De dataset met validatiepunten bevat informatie over de verandering van een
  # cel (0/1) en de aard van de verandering (klasse A -> klasse B). Voor de
  # validatie werden achteraf een aantal moeilijk te onderscheiden klassen
  # samengenomen. De analyse kan dus op 3 niveaus uitgevoerd worden:
  #   (1) verandering - geen verandering, (2) originele landgebruiksklassen en
  # (3) de geaggregeerde klassen. In het onderstaande script worden deze
  # drie validatiesets apart aangemaakt.

  # 1. Puntenset -----
  points <- punten %>%
    na_if("<Null>") %>%
    # Alle <Null> omzetten in NA
    gather(key = klasse, value = oordeel, X2013:change.9) %>%
    # Van breed formaat naar lang formaat
    separate(col = klasse, into = c("klasse", "eval"), sep = "\\.") %>%
    replace_na(list(eval = 0)) %>%
    # Voeg evaluator toe
    drop_na() %>%
    # Alle NA laten vallen
    rename(
      x = POINT_X, y = POINT_Y, lu2013 = change2013, lu2016 = change2016,
      verandering = type
    ) %>%
    mutate(oordeel = recode(oordeel, "Nee" = "nochange", "Ja" = "change")) %>%
    mutate(klasse = gsub("\\..*", "", klasse)) %>%
    # Alles na "." weglaten -> \\.. definieert . en * betekent "alles na"
    mutate(klasse = recode(klasse,
                           X2013 = "lu2013", X2016 = "lu2016",
                           change = "verandering"
    )) %>%
    mutate(lu2013 = recode(lu2013,
                           "1" = "Open natuur", "2" = "Bos", "3" = "Grasland", "4" = "Akker",
                           "5" = "Urbaan",
                           "6" = "Laag groen", "7" = "Hoog groen", "8" = "Water", "9" = "Overig"
    )) %>%
    # codes naar tekst
    mutate(lu2016 = recode(lu2016,
                           "1" = "Open natuur", "2" = "Bos", "3" = "Grasland",
                           "4" = "Akker", "5" = "Urbaan", "6" = "Laag groen",
                           "7" = "Hoog groen", "8" = "Water", "9" = "Overig"
    )) %>%
    # codes naar tekst
    rowwise() %>%
    mutate(lu_c = ifelse(klasse == "lu2013" & lu2013 == oordeel,
                         1, ifelse(klasse == "lu2016" & lu2016 == oordeel,
                                   1, 0
                         )
    )) %>%
    # Check of de gevalideerde landgebruiken overeenkomen met de LG van de kaart
    mutate(change_c = ifelse(verandering == oordeel,
                             1, 0
    )) %>%
    # Check of de beoordeling "change/nochange" overeenkomt met die van de
    # LG-kaart
    ungroup() %>%
    group_by(objectid) %>%
    mutate(n = n() / 3) %>%
    ungroup() %>%
    filter(n >= 1) %>%
    # alleen punten die gevalideerd zijn
    arrange(objectid) %>%
    left_join(dplyr::select(tbltrans, valid, lu), by = c("lu2013" = "lu")) %>%
    rename(luval13 = valid) %>%
    # validatieklassen toevoegen (= aggregatie van oorspronkelijke lu-klassen)
    left_join(dplyr::select(tbltrans, valid, lu), by = c("lu2016" = "lu")) %>%
    rename(luval16 = valid) %>%
    left_join(dplyr::select(tbltrans, valid, lu), by = c("oordeel" = "lu")) %>%
    rename(oordeelval = valid) %>%
    group_by(objectid, eval) %>%
    mutate(oordeelval = ifelse(is.na(oordeelval),
                               ifelse(identical(oordeelval[1], oordeelval[2]),
                                      "nochange", "change"
                               ), oordeelval
    )) %>%
    # Aanpassen beoordeling "verandering" -> als de validatieklasse 2 x hetzelfde
    # is per evaluator, dan "nochange"
    rowwise() %>%
    mutate(veranderingval = ifelse(luval13 == luval16, "nochange", "change")) %>%
    mutate(luval_c = ifelse(klasse == "lu2013" & luval13 == oordeelval, 1,
                            ifelse(klasse == "lu2016" & luval16 == oordeelval, 1, 0)
    )) %>%
    # Check of de gevalideerde landgebruiken overeenkomen met de LG van de kaart
    mutate(changeval_c = ifelse(veranderingval == oordeelval, 1, 0)) %>%
    # Check of de beoordeling "change/nochange" overeenkomt met die van
    # de LG-kaart
    ungroup()

  points_id <- points %>%
    # puntenset met 1 waarde per objectid
    # filter(klasse != "verandering") %>%
    spread(klasse, oordeelval) %>%
    group_by(objectid, eval) %>%
    summarise(
      luval13 = first(luval13), luval16 = first(luval16),
      verand = first(veranderingval),
      lu13oord = first(na.omit(lu2013)), lu16oord = first(na.omit(lu2016)),
      verandoord = first(na.omit(verandering))
    ) %>%
    group_by(objectid) %>%
    sample_n(1) %>% #1 random classificatie w gekozen bij conflict tss experten
    na.omit() %>%
    mutate(codeval13 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", luval13,
                            perl = TRUE
    )) %>%
    mutate(codeval16 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", luval16,
                            perl = TRUE
    )) %>%
    unite(changeclass, codeval13, codeval16, sep = "_") %>%
    mutate(codeval13 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", lu13oord,
                            perl = TRUE
    )) %>%
    mutate(codeval16 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", lu16oord,
                            perl = TRUE
    )) %>%
    unite(changeclassref, codeval13, codeval16, sep = "_") %>%
    left_join(punten[, c('objectid', 'POINT_X', 'POINT_Y')],
              by = c("objectid" = "objectid")) %>%
    mutate(
      changeclass = as.factor(changeclass),
      changeclassref = as.factor(changeclassref)
    ) %>%
    filter(#No water classes
      !changeclass %in% c("W_W", "ON_W", "O_W"),
      !changeclassref %in% c("W_W", "ON_W", "O_W")
    ) %>%
    droplevels() %>%
    left_join(unique(tbltrans[,c("valid", "valid_eng")]),
              by = c("lu13oord" = "valid")) %>%
    rename(valid_eng = lu13oord) %>%
    left_join(unique(tbltrans[,c("valid", "valid_eng")]),
              by = c("lu13oord" = "valid")) %>%
    rename(lu13oord_eng = valid_eng ) %>%
    left_join(unique(tbltrans[,c("valid", "valid_eng")]),
              by = c("lu16oord" = "valid")) %>%
    rename(lu16oord_eng = valid_eng )
  save(points_id, tbltrans, file = "data/validation.Rdata")

  ############################# Get areas #################################
  lgarea <- combine %>%
    dplyr::select(LG2013_ChangeCla, LG2016_ChangeCla, Count) %>%
    rename(
      lu2013 = LG2013_ChangeCla, lu2016 = LG2016_ChangeCla,
      count = Count
    ) %>%
    mutate(lu2013 = as.factor(lu2013),
           lu2016 = as.factor(lu2016)) %>%
    left_join(dplyr::select(tbltrans, valid, lucode),
              by = c("lu2013" = "lucode")
    ) %>%
    rename(val2013 = valid) %>%
    left_join(dplyr::select(tbltrans, valid, lucode),
              by = c("lu2016" = "lucode")
    ) %>%
    rename(val2016 = valid) %>%
    mutate(codeval13 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", val2013,
                            perl = TRUE
    )) %>%
    mutate(codeval16 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", val2016,
                            perl = TRUE
    )) %>%
    unite(class, codeval13, codeval16, sep = "_") %>%
    dplyr::select(class, count) %>%
    group_by(class) %>%
    summarise(count = sum(count)) %>%
    filter(!str_detect(class, "W")) %>%
    # alles met water in weglaten
    droplevels() %>%
    mutate(area = count / sum(count))

  #change-no change
  lgarea_change <- lgarea %>%
    mutate(change = ifelse(!class %in% c("A_A", "HG_HG", "O_O", "ON_ON", "U_U"),
                           "change", "nochange"
    )) %>%
    group_by(change) %>%
    summarise(count = sum(count), area = sum(area))

  # 2.3. Originele landgebruiksklassen ----

  lgarea_lu <- combine %>%
    dplyr::select(LG2013_ChangeCla, LG2016_ChangeCla, Count) %>%
    rename(
      lu2013code = LG2013_ChangeCla, lu2016code = LG2016_ChangeCla,
      count = Count
    ) %>%
    mutate(lu2013code = as.factor(lu2013code),
           lu2016code = as.factor(lu2016code)) %>%
    left_join(dplyr::select(tbltrans, lucode, lu),
              by = c("lu2013code" = "lucode")
    ) %>%
    rename(lu2013 = lu) %>%
    left_join(dplyr::select(tbltrans, lucode, lu),
              by = c("lu2016code" = "lucode")
    ) %>%
    rename(lu2016 = lu) %>%
    mutate(code13 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", lu2013,
                         perl = TRUE
    )) %>%
    mutate(code16 = gsub("\\b(\\pL)\\pL{2,}|.", "\\U\\1", lu2016,
                         perl = TRUE
    )) %>%
    unite(class, code13, code16, sep = "_") %>%
    dplyr::select(class, count) %>%
    group_by(class) %>%
    summarise(count = sum(count)) %>%
    filter(!str_detect(class, "W")) %>%
    # alles met water in weglaten
    droplevels() %>%
    mutate(area = count / sum(count))
}


