##---------------------------------------------------------------------------
## Author:      Simone Conti
## Copyright:   (C) 2014 Simone Conti 
## License:     GPLv3+
##---------------------------------------------------------------------------



# Macro to create a new wordom test

macro(wrd_add_test wrdtest)
    execute_process( COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/tests/${wrdtest})
    add_test ( 
        NAME ${wrdtest}
        COMMAND ${CMAKE_COMMAND} 
                    -D WORDOM=${CMAKE_BINARY_DIR}/wordom
                    -D SRCDIR=${CMAKE_CURRENT_SOURCE_DIR}/${wrdtest}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/${wrdtest}/runtest.cmake
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/tests/${wrdtest}
    )
endmacro(wrd_add_test)



# Macro to compare the output generated and the ref one in srcdir

macro(wrd_compare_file FNAME)

    message( "Comparing ${SRCDIR}/${FNAME}.ref with ${FNAME}")
    execute_process( COMMAND ${CMAKE_COMMAND} -E compare_files ${SRCDIR}/${FNAME}.ref ${FNAME} RESULT_VARIABLE TEST_RESULT )

    if (TEST_RESULT)
        message(FATAL_ERROR "Failed: The output of ${FNAME} did not match ${SRCDIR}/${FNAME}.ref")
    endif(TEST_RESULT)

endmacro(wrd_compare_file)



# Macro to run wordom with a specific command line args. Redirect out and err in LOGFILE.log and LOGFILE.err

macro(wrd_run WORDOM LOGFILE ARGS)

    message( "Running ${WORDOM} ${ARGS}")

    # Run wordom
    execute_process(
        COMMAND ${WORDOM} ${ARGS}
        OUTPUT_FILE     ${LOGFILE}.log
        ERROR_FILE      ${LOGFILE}.err
        ERROR_VARIABLE  TEST_ERROR
        RESULT_VARIABLE TEST_RESULT
    )

    # Check the exit code of wordom
    if( TEST_RESULT )
        message( FATAL_ERROR "Wordom exited with code ${TEST_ERRROR} != 0" )
    endif( TEST_RESULT )

endmacro(wrd_run)

