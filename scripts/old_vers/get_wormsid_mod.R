function (sci_com, searchtype = "scientific", marine_only = TRUE, 
          fuzzy = NULL, accepted = FALSE, ask = TRUE, messages = TRUE, 
          rows = NA, query = NULL, ...) 
{
  assert(sci_com, c("character", "taxon_state"))
  assert(searchtype, "character")
  assert(marine_only, "logical")
  assert(fuzzy, "logical")
  assert(accepted, "logical")
  assert(ask, "logical")
  assert(messages, "logical")
  assert_rows(rows)
  pchk(query, "sci_com")
  if (inherits(sci_com, "character")) {
    tstate <- taxon_state$new(class = "wormsid", names = sci_com)
    items <- sci_com
  }
  else {
    assert_state(sci_com, "wormsid")
    tstate <- sci_com
    sci_com <- tstate$taxa_remaining()
    items <- c(sci_com, tstate$taxa_completed())
  }
  prog <- progressor$new(items = items, suppress = !messages)
  done <- tstate$get()
  for (i in seq_along(done)) prog$completed(names(done)[i], 
                                            done[[i]]$att)
  prog$prog_start()
  for (i in seq_along(sci_com)) {
    direct <- FALSE
    mssg(messages, "\nRetrieving data for taxon '", sci_com[i], 
         "'\n")
    if (!searchtype %in% c("scientific", "common")) {
      stop("'searchtype' must be one of 'scientific' or 'common'", 
           call. = FALSE)
    }
    wmdf <- switch(searchtype, scientific = worms_worker(sci_com[i], 
                                                         worrms::wm_records_name, rows, marine_only, fuzzy %||% 
                                                           TRUE, ...), common = worms_worker(sci_com[i], 
                                                                                             worrms::wm_records_common, rows, marine_only, fuzzy %||% 
                                                                                               FALSE, ...))
    mm <- NROW(wmdf) > 1
    if (!inherits(wmdf, "tbl_df") || NROW(wmdf) == 0) {
      wmid <- NA_character_
      att <- "not found"
      wmdf <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                NA, NA, NA, NA, NA, NA, NA, NA, NA)
    }
    else {
      wmdf <- suppressWarnings(data.frame(wmdf))
      wmdf <- wmdf[, c("AphiaID", "rank", "status", "unacceptreason", 
                       "valid_AphiaID", "valid_name", "kingdom", "phylum", 
                       "class", "order", "family", "genus", "isMarine", 
                       "isBrackish", "isFreshwater", "isTerrestrial", 
                       "isExtinct")]
      names(wmdf)[1] <- "id"
      if (accepted) {
        wmdf <- wmdf[wmdf$status %in% "accepted", ]
      }
      wmdf <- sub_rows(wmdf, rows)
      if (nrow(wmdf) == 0) {
        mssg(messages, m_not_found_sp_altclass)
        wmid <- NA_character_
        att <- "not found"
        wmdf <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, 
                  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
      }
      if (nrow(wmdf) == 1) {
        wmid <- wmdf$id
        att <- "found"
      }
      if (nrow(wmdf) > 1) {
        names(wmdf)[grep("scientificname", names(wmdf))] <- "target"
        matchtmp <- wmdf[tolower(wmdf$target) %in% tolower(sci_com[i]), 
                         "id"]
        if (length(matchtmp) == 1) {
          wmid <- matchtmp
          direct <- TRUE
          att <- "found"
        }
        else {
          wmid <- NA_character_
          att <- "not found"
        }
      }
      if (any(nrow(wmdf) > 1 && is.na(wmid) | nrow(wmdf) > 
              1 && att == "found" & length(wmid) > 1)) {
        if (ask) {
          names(wmdf)[grep("scientificname", names(wmdf))] <- "target"
          wmdf <- wmdf[order(wmdf$target), ]
          message("\n\n")
          print(wmdf)
          message("\nMore than one WORMS ID found for taxon '", 
                  sci_com[i], "'!\n\n                  Enter rownumber of taxon (other inputs will return 'NA'):\n")
          take <- scan(n = 1, quiet = TRUE, what = "raw")
          if (length(take) == 0) {
            take <- "notake"
            att <- "nothing chosen"
          }
          if (take %in% seq_len(nrow(wmdf))) {
            take <- as.numeric(take)
            message("Input accepted, took taxon '", as.character(wmdf$target[take]), 
                    "'.\n")
            wmid <- wmdf$id[take]
            att <- "found"
          }
          else {
            wmid <- NA_character_
            mssg(messages, "\nReturned 'NA'!\n\n")
            att <- "not found"
          }
        }
        else {
          if (length(wmid) != 1) {
            warning(sprintf(m_more_than_one_found, "Worms ID", 
                            sci_com[i]), call. = FALSE)
            wmid <- NA_character_
            att <- m_na_ask_false
          }
        }
      }
    }
    res <- list(id = as.character(wmid), att = att, multiple = mm, 
                direct = direct, rank = as.character(wmdf[2]), status = as.character(wmdf[3]), 
                unaccreason = as.character(wmdf[4]), validid = as.character(wmdf[5]), 
                validname = as.character(wmdf[6]), kingdom = as.character(wmdf[7]), 
                phylum = as.character(wmdf[8]), class = as.character(wmdf[9]), 
                order = as.character(wmdf[10]), family = as.character(wmdf[11]), 
                genus = as.character(wmdf[12]), marine = as.character(wmdf[13]), 
                brackish = as.character(wmdf[14]), fresh = as.character(wmdf[15]), 
                terr = as.character(wmdf[16]), dead = as.character(wmdf[17]))
    prog$completed(sci_com[i], att)
    prog$prog(att)
    tstate$add(sci_com[i], res)
  }
  out <- tstate$get()
  ids <- data.frame(query = sci_com, aphia_id = pluck_un(out, 
                                                         "id", ""), rank = pluck_un(out, "rank", ""), status = pluck_un(out, 
                                                                                                                        "status", ""), unacceptreason = pluck_un(out, "unaccreason", 
                                                                                                                                                                 ""), valid_aphiaid = pluck_un(out, "validid", ""), valid_name = pluck_un(out, 
                                                                                                                                                                                                                                          "validname", ""), kingdom = pluck_un(out, "kingdom", 
                                                                                                                                                                                                                                                                               ""), phylum = pluck_un(out, "phylum", ""), class = pluck_un(out, 
                                                                                                                                                                                                                                                                                                                                           "class", ""), order = pluck_un(out, "order", ""), family = pluck_un(out, 
                                                                                                                                                                                                                                                                                                                                                                                                               "family", ""), genus = pluck_un(out, "genus", ""), ismarine = pluck_un(out, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      "marine", ""), isbrackish = pluck_un(out, "brackish", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           ""), isfreshwater = pluck_un(out, "fresh", ""), isterrestrial = pluck_un(out, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    "terr", ""), isextinct = pluck_un(out, "dead", ""), found = pluck_un(out, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         "att", ""), pattern_match = pluck_un(out, "direct", logical(1)))
  on.exit(prog$prog_summary(), add = TRUE)
  on.exit(tstate$exit, add = TRUE)
  add_uri(ids, get_url_templates$worms)
}
