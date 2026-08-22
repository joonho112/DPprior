% Please edit documentation in R/01_stirling.R
\name{.stirling_cache_info}
\alias{.stirling_cache_info}
\title{Inspect the Session-local Stirling Cache}
\usage{
.stirling_cache_info()
}
\value{
A list containing the cache policy, largest cached J, hit/miss
counts, and retained bytes.
}
\description{
Internal diagnostic helper used to audit cache provenance and memory.
}
\keyword{internal}
