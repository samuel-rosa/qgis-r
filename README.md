---
output: html_document
---

## May 2017

To avoid errors in QGIS-R integration due to 'user/local/lib/R/site-library is not writable', 
edit `etc/R/Renviron`.
http://stackoverflow.com/questions/2698269/how-do-you-change-library-location-in-r
https://cran.r-project.org/doc/manuals/r-devel/R-admin.html#Managing-libraries
http://stackoverflow.com/questions/15170399/changing-r-default-library-path-using-libpaths-in-rprofile-site-fails-to-work
    
    R_LIBS_SITE=${R_LIBS_SITE-'/home/alessandro/Rlibs:/usr/local/lib/R/site-library:/usr/lib/R/site-library:/usr/lib/R/library'}

## More info at
https://docs.qgis.org/2.18/en/docs/user_manual/processing/scripts.html

https://docs.qgis.org/2.18/en/docs/user_manual/processing/3rdParty.html#r-creating-r-scripts

https://docs.qgis.org/2.18/en/docs/training_manual/processing/r_intro.html

https://docs.qgis.org/2.18/en/docs/training_manual/processing/r_syntax.html#r-syntax

# Available Learners:
https://topepo.github.io/caret/available-models.html
