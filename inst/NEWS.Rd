\name{frmanews}
\title{fRMA News}
\encoding{UTF-8}

\section{Version 1.19}{
  \itemize{
    \item Correct citation for use of frma with Gene and Exon ST arrays
    added. 
    \item The internal workings of the frma function have been modified
    to ease the eventual transition AffyBatch (affy) objects to FeatureSet
    (oligo) objects as input. 
    \item Changed the way in which the user can pass their own fRMA
    vectors to the frma function via the "input.vecs"
    argument. Previously, the user could supply a partial set of vectors
    and the remaining ones would be obtained from the default frmavecs
    packages. In retrospect, this seems very likely to introduce a
    mistake. Now, the user has to either supply all of their own vectors
    or use all of the default vectors. If anyone was mixing vectors
    before, this can still be accomplished by loading the respective
    frmavecs package and manually creating your own mixture of vectors. 
    \item Bug fix: on certain operating systems, due to differences in
    the locale, the probe ordering returned by the pm function (affy) did
    not match the order stored in the frmavecs packages. Updated frmavecs
    packages will now contain the pmindex values for all probes and affy
    ids for all probesets. Also the frma and barcode functions now check
    that the order matches. If the order is not the same, these functions
    attempt to reorder the probes to match.   
    \item Added NEWS.Rd file.
  }
}
