import_cl_results <- function (cl_file) {

    ### read file
    conn <- file(cl_file, open='r')
    lines <- readLines(conn)
    close(conn)


    ### build table
    DF <- data.frame(NULL, stringsAsFactors = FALSE)

    for (x in lines) {

        if (grepl('^#', x)) { next }

        v <- strsplit(x, '\\t')[[1]]

        if (grepl('^>Cluster', x)) {
    
            cl_repres <- v[2]
            id <- v[2]
            len <- v[3]

        } else {

            id <- v[1]
            len <- v[5]

        }

        df <- data.frame(cl_repres, cl_member = id, cl_member_len = len, stringsAsFactors = FALSE)
        DF <- rbind(DF, df)
    }


    ### return table
    return(DF)

}
