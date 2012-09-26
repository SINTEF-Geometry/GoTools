IF(WIN32)
#    SET(CTEST_CUSTOM_POST_TEST "type Testing/Temporary/LastTest.log")
    SET(CTEST_CUSTOM_POST_TEST "cd")
ENDIF(WIN32)
