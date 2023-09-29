# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project should adhere to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

### Changed

### Fixed

- Bug fix in `createISO31662Grid` when indexing natural earth data using iso_3166_2 column

### Removed

## [2.0.0] - 2021-10-01

### Added

- New default RUN script of `RUN_2_geospatial_uncor`
- `CreateTradespacePairingGeo` is a helper function crete a tradespace of different combinations of input variables to `createEncounters_2`
- `trackTimetable2NEU` converts tracks in timetable format to an array corresponding to time, north, east, and up coordinates
- `adjustTrack` will adjust tracks based on output of `samplespeedalt`
- `estimateSpeedFromGeodetic` is a helper function to estimate speed based on latitude and longitude
- `preallocAnchors`, `createISO31662Grid`, and `parseGeoTrajDirectory` helper functions for step one
- Persistent system environment variable of `AEM_DIR_GEOPAIR` established

### Changed

- Inputs to `findPairs_1` changed to take advantage of calling `preallocAnchors` within `findPairs_1`. The step one RUN script was updated accordingly.
- `createEncounters_2` no longer reads in pre-generated files of tracks generated from the Bayesian encounter models; instead encounter model based tracks are generated dynamically using the new `UncorEncounterModel` class added to `em-model-manned-bayes`
- `createEncounters_2` will try multiple times to create an encounter with aircraft 1 and will now resample an encounter model to generate new Bayes tracks if needed. The previous version gave up too easily.
- `createEncounters_2` no longer organized loops by clusters and all clustering has been removed
- `createEncounters_2` loads FAA digital obstacle file using `gridDOF` from `em-core`
- `createEncounters_2` stores tracks as a `timetable`, a type of table that associates a time with each row. Helper functions like `findconflict`, `samplespeedalt`, and `findCropIdx` have been updated accordingly
- `createEncounters_2` outputs a table containing encounter metadata
- `createEncounters_2` and `loadTrack` will enforce the altitude range (altRange1_ft_agl, altRange2_ft_agl) when sampling the Bayes tracks. For Bayes tracks, altitude still cannot resampled using `samplespeedalt` and `adjustTrack`.
- `findconflict` now considers HMD, VMD, and time criteria. It previously only assessed HMD.
- `findCropIdx` also now checks for initial horizontal and vertical separation criteria
- Renamed `neu_to_wp_struct` to to `neu2Waypoints`
- Renamed `startup` to `startup_pairing_geo`
- DEM handling updated based on improvements to `msl2agl` from `em-core`
- Better random seed control
- Improved plotting
- Improved organization through the use of the code directory

### Fixed

- `samplespeedalt` handles rare case where the altitude span is less than 25 resulting in `adjustSpan_ft_agl` being empty
- Fixed bug in `findConflict` when calculating `hmd_ft` and `vmd_ft` but tracks are not near each other
- Various functions renamed using camelCase

### Removed

- `croptracks`, `calcLegsTime`, `neu2wpstruct`, `updatewaypointstruct` removed due to the addition of `trackTimetable2NEU` and updates to `findconflict` and `findCropIdx`
- `interpTime` removed because time based interpolation is handled using the built-in functionality of the `timetable` type
- `RUN_2_sUAS_uncor` and `RUN_2_sUAS_unconv` removed because they were deprecated by the more generalized `RUN_2_geospatial_uncor` script
- `RUN_2_sUAS_HAA` removed due to change of sampling Bayesian networks at runtime. The Bayesian helicopter air ambulance model hasn't been tested yet with the `EncounterModel` objects in em-model-manned-bayes. This functionality maybe reintroduced in a future release.

## [1.1.0] - 2020-08-10

### Added

- SPDX headers

### Changed

- Improving handling of Bayes tracks from em-model-manned-bayes

### Fixed

- Addressed bug `placeTrack` where tracks were not correctly rotated

## [1.0.0] - 2020-06-24

### Added

- Initial public release
