#include <check.h>

#include "functions.c"

START_TEST(test_Rsinh)
{
    ck_assert(abs(Rsinh(5)-3.6286041)<=3.6286041*FLT_EPSILON);
    ck_assert(abs(Rsinh(8)-1490.47882579)<=1490.47882579*FLT_EPSILON);

}
END_TEST

START_TEST(test_Rlog)
{
    ck_assert(abs(Rlog(5)-1.60943791)<=1.60943791*FLT_EPSILON);
    ck_assert(abs(Rlog(8)-2.07944154)<=2.07944154*FLT_EPSILON);

}
END_TEST

START_TEST(test_ocp_anode)
{
    ck_assert(abs(ocp_anode(500,800)-4.51681399)<=4.51681399*FLT_EPSILON);
    ck_assert(abs(ocp_anode(800,1000)-4.51681399)<=4.51681399*FLT_EPSILON);

}
END_TEST

START_TEST(test_ocp_cathode)
{
    ck_assert(abs(ocp_cathode(500,800)-0.21864040)<=0.21864040*FLT_EPSILON);
    ck_assert(abs(ocp_cathode(700,1300)-0.21864040)<=0.21864040*FLT_EPSILON);

}
END_TEST

int main(void)
{
    Suite *s1 = suite_create("Core");
    TCase *tc1_1 = tcase_create("Core");
    SRunner *sr = srunner_create(s1);
    int nf;

    
    tcase_add_test(tc1_1, test_Rsinh);
    tcase_add_test(tc1_1, test_Rlog);
    tcase_add_test(tc1_1, test_ocp_anode);
    tcase_add_test(tc1_1, test_ocp_cathode);
    suite_add_tcase(s1, tc1_1);
    
    
    srunner_run_all(sr, CK_ENV);
    nf = srunner_ntests_failed(sr);
    srunner_free(sr);

    return nf == 0 ? 0 : 1;
}
