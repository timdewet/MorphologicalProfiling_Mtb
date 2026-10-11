# =============================================================================
# CLARIOstar plate-reader helpers shared by the Plate60 scripts.
#   read_clariostar_od()  - OD600 per well straight from a .RUC file
#   find_ruc()            - the .RUC for a timepoint ("day 4", "Day7", ...)
#   add_plate_layout()    - reader well -> imaging well, blank-subtracted OD
# =============================================================================

# CLARIOstar .RUC = zip holding an "Absolute Database" file. Each well-channel
# record contains: 0xE0, well label ("A01".."H12"), ..., wavelength (double,
# +18 bytes) and reading (double, +26 bytes). Channel 600 nm is the OD; the
# 900/993 nm channels are pathlength-correction reads and are ignored.
read_clariostar_od <- function(ruc, wavelength = 600) {
  tmp <- tempfile(); on.exit(unlink(tmp, recursive = TRUE))
  utils::unzip(ruc, exdir = tmp)
  f <- list.files(tmp, pattern = "\\.abs$", full.names = TRUE)
  b <- readBin(f, "raw", file.info(f)$size)
  is_row <- b >= as.raw(0x41) & b <= as.raw(0x48)          # "A".."H"
  is_dig <- b >= as.raw(0x30) & b <= as.raw(0x39)
  n <- length(b)
  hit <- which(b[1:(n - 34)] == as.raw(0xE0))
  hit <- hit[is_row[hit + 1] & is_dig[hit + 2] & is_dig[hit + 3] & b[hit + 4] == as.raw(0)]
  rd  <- function(at) vapply(at, function(i) readBin(b[i:(i + 7)], "double", 1, 8, endian = "little"), 0)
  out <- data.table(well  = vapply(hit, function(i) rawToChar(b[(i + 1):(i + 3)]), ""),
                    wl    = rd(hit + 18),
                    value = rd(hit + 26))
  out <- unique(out[wl == wavelength, .(well, value)])
  stopifnot(nrow(out) == 96, !anyDuplicated(out$well))
  out[, `:=`(row = substr(well, 1, 1), col = as.integer(substring(well, 2)))]
  out[]
}

# .RUC files in data_dir whose name carries the timepoint ("day 4", "Day7";
# case-insensitive).  Returns all matches; callers decide what to do with >1.
find_ruc <- function(data_dir, tag) {
  day <- sub("^day", "", regmatches(tag, regexpr("day\\d+", tag)))
  cand <- list.files(data_dir, pattern = "\\.RUC$", full.names = TRUE)
  if (length(day)) cand[grepl(paste0("day ?", day, "_"), basename(cand), ignore.case = TRUE)] else cand
}

# Reader block A2-F11 = imaging wells A1-F10 with columns mirrored (reader A2
# -> imaging A10, reader A11 -> imaging A1).  Row H stays flat (~0.31) on every
# day and is used as the blank: the other edge wells (column 1, column 12,
# row G) grow over the week, so they are cultures, not blanks.
BLANK_ROW <- "H"
add_plate_layout <- function(od) {
  od[, img_well := fifelse(row %in% LETTERS[1:6] & col %in% 2:11, paste0(row, 12L - col), NA_character_)]
  od[, blank := median(value[row == BLANK_ROW])]
  od[, od_net := value - blank]
  od[]
}
