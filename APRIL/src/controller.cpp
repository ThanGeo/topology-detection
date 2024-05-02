#include "controller.h"

namespace APRIL
{
    spatial_lib::QueryT* g_query;

    void setupAPRILIntermediateFilter(spatial_lib::QueryT *query) {
        // store query info
        g_query = query;
    }
    
    namespace optimized
    {
        namespace relate
        {
            /**
             * WARNING: standard and optimized relate functions are THE SAME right now, but I keep separated because future optimization might happen.
             * 
            */

            static void APRILIntersectionFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = intersectionAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineIntersectionJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILInsideFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = insideAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineInsideJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILDisjointFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = disjointAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineDisjointJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILEqualFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = equalAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineEqualJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILMeetFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = meetAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                
                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineMeetJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILContainsFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = containsAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineContainsJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILCoversFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = coversAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                
                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineCoversJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILCoveredByFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = coveredByAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineCoveredByJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }
        }


        clock_t temp_timer;
        namespace find_relation
        {
            static void specializedTopologyRinSContainment(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                        iFilterResult = RinSContainmentAPRILUncompressed(aprilR, aprilS);
                    } else {
                        // compressed
                        iFilterResult = compressed::RinSContainmentAPRIL(aprilR, aprilS);
                    }
                    
                    
                    // switch based on result
                    switch(iFilterResult) {
                        // true hits, count and return
                        case spatial_lib::TR_INSIDE:
                        case spatial_lib::TR_DISJOINT:
                        case spatial_lib::TR_INTERSECT:
                            // result
                            spatial_lib::countTopologyRelationResult(iFilterResult);

                            // if(iFilterResult == spatial_lib::TR_INTERSECT){
                            //     // printf("%u,%u\n", idR, idS);
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     if(vertexCounts.first > 2000 && vertexCounts.second > 2000) {
                            //         printf("polygons %u,%u vertex counts %u,%u\n", idR, idS, vertexCounts.first, vertexCounts.second);
                            //     }
                            // }

                            // if(iFilterResult == spatial_lib::TR_INTERSECT){
                            //     printf("%u,%u\n", idR, idS);
                            // }

                            // if (idR == 220594 && idS == 357429) {
                            // if (idR == 206434 && idS == 9173759) {
                            //     printf("Time passed evaluating pair %u,%u: %f seconds\n", idR, idS, spatial_lib::time::getElapsedTime(temp_timer));
                            //     printf("A/F intervals: %u/%u and %u/%u\n", aprilR->numIntervalsALL, aprilR->numIntervalsFULL, aprilS->numIntervalsALL, aprilS->numIntervalsFULL);
                            // }


                            // dev-pol-stats - true hit
                            // auto itR = spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.find(idR);
                            // if (itR == spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.end()) {
                            //     // count as first-time TH
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.insert(std::make_pair(idR, 1));
                            //     // count vertices
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapR.insert(std::make_pair(idR, vertexCounts.first));
                            // } else {
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapR[idR]++;
                            // }
                            // auto itS = spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.find(idS);
                            // if (itS == spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.end()) {
                            //     // count as first-time TH
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.insert(std::make_pair(idS, 1));
                            //     // count vertices
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapS.insert(std::make_pair(idS, vertexCounts.second));
                            // } else {
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapS[idS]++;
                            // }




                            // time
                            spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                            return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                int relation;
                // switch based on result
                switch(iFilterResult) {            
                    // inconclusive, do selective refinement
                    // result holds what type of refinement needs to be made
                    case spatial_lib::REFINE_INSIDE_COVEREDBY_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineInsideCoveredbyTruehitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_DISJOINT_INSIDE_COVEREDBY_MEET_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineDisjointInsideCoveredbyMeetIntersect(idR, idS);
                        break;
                    
                }

                // dev-pol-stats - refinement
                // auto itR = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.find(idR);
                // if (itR == spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.end()) {
                //     // count as first-time RF
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.insert(std::make_pair(idR, 1));
                //     // count vertices
                //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapR.insert(std::make_pair(idR, vertexCounts.first));
                // } else {
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR[idR]++;
                // }
                // auto itS = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.find(idS);
                // if (itS == spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.end()) {
                //     // count as first-time RF
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.insert(std::make_pair(idS, 1));
                //     // count vertices
                //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapS.insert(std::make_pair(idS, vertexCounts.second));
                // } else {
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS[idS]++;
                // }


                // time
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
                // count the result
                spatial_lib::countTopologyRelationResult(relation);
                // if(relation == spatial_lib::TR_INTERSECT){
                //     printf("%u,%u\n", idR, idS);
                // }
            }

            static void specializedTopologySinRContainment(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                        iFilterResult = SinRContainmentAPRILUncompressed(aprilR, aprilS);
                    } else {
                        // compressed
                        iFilterResult = compressed::SinRContainmentAPRIL(aprilR, aprilS);
                    }
                    
                    // switch based on result
                    switch(iFilterResult) {
                        // true hits, count and return
                        case spatial_lib::TR_CONTAINS:
                        case spatial_lib::TR_DISJOINT:
                        case spatial_lib::TR_INTERSECT:
                            // if(iFilterResult == spatial_lib::TR_INTERSECT){
                            //     printf("%u,%u\n", idR, idS);
                            // }






                            // dev-pol-stats - true hit
                            // auto itR = spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.find(idR);
                            // if (itR == spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.end()) {
                            //     // count as first-time TH
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.insert(std::make_pair(idR, 1));
                            //     // count vertices
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapR.insert(std::make_pair(idR, vertexCounts.first));
                            // } else {
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapR[idR]++;
                            // }
                            // auto itS = spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.find(idS);
                            // if (itS == spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.end()) {
                            //     // count as first-time TH
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.insert(std::make_pair(idS, 1));
                            //     // count vertices
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapS.insert(std::make_pair(idS, vertexCounts.second));
                            // } else {
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapS[idS]++;
                            // }






                            // result
                            spatial_lib::countTopologyRelationResult(iFilterResult);
                            // time
                            spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                            return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                int relation;
                // switch based on result
                switch(iFilterResult) {            
                    // inconclusive, do selective refinement
                    // result holds what type of refinement needs to be made
                    case spatial_lib::REFINE_CONTAINS_COVERS_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineContainsCoversTruehitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_DISJOINT_CONTAINS_COVERS_MEET_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineDisjointContainsCoversMeetIntersect(idR, idS);
                        break;
                    
                }


                // dev-pol-stats - refinement
                // auto itR = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.find(idR);
                // if (itR == spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.end()) {
                //     // count as first-time RF
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.insert(std::make_pair(idR, 1));
                //     // count vertices
                //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapR.insert(std::make_pair(idR, vertexCounts.first));
                // } else {
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR[idR]++;
                // }
                // auto itS = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.find(idS);
                // if (itS == spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.end()) {
                //     // count as first-time RF
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.insert(std::make_pair(idS, 1));
                //     // count vertices
                //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapS.insert(std::make_pair(idS, vertexCounts.second));
                // } else {
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS[idS]++;
                // }






                // time
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
                // count the result
                spatial_lib::countTopologyRelationResult(relation);
                // if(relation == spatial_lib::TR_INTERSECT){
                //     printf("%u,%u\n", idR, idS);
                // }
            }

            static void specializedTopologyIntersection(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                        iFilterResult = MBRIntersectionAPRILUncompressed(aprilR, aprilS);
                    } else {
                        // compressed
                        iFilterResult = compressed::MBRIntersectionAPRIL(aprilR, aprilS);
                    }
                    
                   
                   
                   
                    // switch based on result
                    switch(iFilterResult) {
                        // true hits, count and return
                        case spatial_lib::TR_DISJOINT:
                        case spatial_lib::TR_INTERSECT:
                            // if(iFilterResult == spatial_lib::TR_INTERSECT){
                            //     printf("%u,%u\n", idR, idS);
                            // }

                            // dev-pol-stats - true hit
                            // auto itR = spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.find(idR);
                            // if (itR == spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.end()) {
                            //     // count as first-time TH
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.insert(std::make_pair(idR, 1));
                            //     // count vertices
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapR.insert(std::make_pair(idR, vertexCounts.first));
                            // } else {
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapR[idR]++;
                            // }
                            // auto itS = spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.find(idS);
                            // if (itS == spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.end()) {
                            //     // count as first-time TH
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.insert(std::make_pair(idS, 1));
                            //     // count vertices
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapS.insert(std::make_pair(idS, vertexCounts.second));
                            // } else {
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapS[idS]++;
                            // }






                            // result
                            spatial_lib::countTopologyRelationResult(iFilterResult);
                            // time
                            spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                            return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);


                // dev-pol-stats - refinement
                // auto itR = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.find(idR);
                // if (itR == spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.end()) {
                //     // count as first-time RF
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.insert(std::make_pair(idR, 1));
                //     // count vertices
                //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapR.insert(std::make_pair(idR, vertexCounts.first));
                // } else {
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR[idR]++;
                // }
                // auto itS = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.find(idS);
                // if (itS == spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.end()) {
                //     // count as first-time RF
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.insert(std::make_pair(idS, 1));
                //     // count vertices
                //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapS.insert(std::make_pair(idS, vertexCounts.second));
                // } else {
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS[idS]++;
                // }





                // time 
                spatial_lib::time::markRefinementFilterTimer();
                // refine
                int relation = spatial_lib::refineDisjointMeetIntersect(idR, idS);
                // time
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
                // count the result
                spatial_lib::countTopologyRelationResult(relation);
                // if(relation == spatial_lib::TR_INTERSECT){
                //     printf("%u,%u\n", idR, idS);
                // }
            }

            static void specializedTopologyEqual(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                        iFilterResult = equalMBRsAPRILUncompressed(idR, idS, aprilR, aprilS);
                    } else {
                        iFilterResult = compressed::equalMBRsAPRIL(idR, idS, aprilR, aprilS);
                    }

                    // if(idR == 4320271 && idS == 5510423) {
                    //     printf("%u,%u ifilter result: %d\n", idR, idS, iFilterResult);
                    //     exit(0);
                    // }
                    
                
                    // switch based on result
                    switch(iFilterResult) {
                        // true hits, count and return
                        case spatial_lib::TR_DISJOINT:
                        // case spatial_lib::TR_EQUAL:
                        case spatial_lib::TR_COVERED_BY:
                        case spatial_lib::TR_COVERS:
                        case spatial_lib::TR_INTERSECT:
                            // if(iFilterResult == spatial_lib::TR_INTERSECT) {
                            //     printf("%u,%u\n",idR,idS);
                            // }


                            // dev-pol-stats - true hit
                            // auto itR = spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.find(idR);
                            // if (itR == spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.end()) {
                            //     // count as first-time TH
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapR.insert(std::make_pair(idR, 1));
                            //     // count vertices
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapR.insert(std::make_pair(idR, vertexCounts.first));
                            // } else {
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapR[idR]++;
                            // }
                            // auto itS = spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.find(idS);
                            // if (itS == spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.end()) {
                            //     // count as first-time TH
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapS.insert(std::make_pair(idS, 1));
                            //     // count vertices
                            //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                            //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapS.insert(std::make_pair(idS, vertexCounts.second));
                            // } else {
                            //     spatial_lib::g_queryOutput.polygonStats.THappearancesMapS[idS]++;
                            // }







                            // result
                            spatial_lib::countTopologyRelationResult(iFilterResult);
                            // time
                            spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                            return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                int relation;
                // switch based on result
                switch(iFilterResult) {            
                    // inconclusive, do selective refinement
                    // result holds what type of refinement needs to be made
                    case spatial_lib::REFINE_COVEREDBY_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineCoveredbyTrueHitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_COVERS_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineCoversTrueHitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_COVERS_COVEREDBY_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineCoversCoveredByTrueHitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_EQUAL_COVERS_COVEREDBY_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineEqualCoversCoveredByTrueHitIntersect(idR, idS);
                        break;
                }

                // dev-pol-stats - refinement
                // auto itR = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.find(idR);
                // if (itR == spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.end()) {
                //     // count as first-time RF
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR.insert(std::make_pair(idR, 1));
                //     // count vertices
                //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapR.insert(std::make_pair(idR, vertexCounts.first));
                // } else {
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapR[idR]++;
                // }
                // auto itS = spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.find(idS);
                // if (itS == spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.end()) {
                //     // count as first-time RF
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS.insert(std::make_pair(idS, 1));
                //     // count vertices
                //     std::pair<uint,uint> vertexCounts = spatial_lib::getVertexCountsOfPair(idR,idS);
                //     spatial_lib::g_queryOutput.polygonStats.polVerticesMapS.insert(std::make_pair(idS, vertexCounts.second));
                // } else {
                //     spatial_lib::g_queryOutput.polygonStats.RFappearancesMapS[idS]++;
                // }








                // time
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
                // count the result
                spatial_lib::countTopologyRelationResult(relation);
                // if(relation == spatial_lib::TR_EQUAL) {
                //     printf("%u,%u\n",idR,idS);
                // }
            }

            static void APRILFindRelation(uint idR, uint idS, int relationCase) {
                
                // use appropriate query function
                switch (relationCase) {
                    case spatial_lib::MBR_R_IN_S:
                        specializedTopologyRinSContainment(idR, idS);
                        break;
                    case spatial_lib::MBR_S_IN_R:
                        specializedTopologySinRContainment(idR, idS);
                        break;
                    case spatial_lib::MBR_EQUAL:
                        specializedTopologyEqual(idR, idS);
                        break;
                    case spatial_lib::MBR_INTERSECT:
                        specializedTopologyIntersection(idR, idS);
                        break;
                    default:
                        // not supported/unknown
                        exit(-1);
                        break;
                }
            }
        }

        void IntermediateFilterEntrypoint(uint idR, uint idS, int relationCase) {
            // time 
            spatial_lib::time::markiFilterTimer();
            // count post mbr candidate
            // spatial_lib::g_queryOutput.postMBRFilterCandidates += 1;
            
            // if (idR == 220594 && idS == 357429) {
            // if (idR == 206434 && idS == 9173759) {
            //     temp_timer = spatial_lib::time::getNewTimer();
            // }

            switch (g_query->type) {
                case spatial_lib::Q_INTERSECT:
                    relate::APRILIntersectionFilter(idR, idS);
                    break;
                case spatial_lib::Q_INSIDE:
                    relate::APRILInsideFilter(idR, idS);
                    break;
                case spatial_lib::Q_DISJOINT:
                    relate::APRILDisjointFilter(idR, idS);
                    break;
                case spatial_lib::Q_EQUAL:
                    relate::APRILEqualFilter(idR, idS);
                    break;
                case spatial_lib::Q_MEET:
                    relate::APRILMeetFilter(idR, idS);
                    break;
                case spatial_lib::Q_CONTAINS:
                    relate::APRILContainsFilter(idR, idS);
                    break;
                case spatial_lib::Q_COVERS:
                    relate::APRILCoversFilter(idR, idS);
                    break;
                case spatial_lib::Q_COVERED_BY:   
                    relate::APRILCoveredByFilter(idR, idS);
                    break;
                case spatial_lib::Q_FIND_RELATION:
                    find_relation::APRILFindRelation(idR, idS, relationCase);
                    break;
                default:
                    fprintf(stderr, "Invalid query and APRIL filter combination.\n");
                    exit(-1);
                    break;
            }
        }   
    }





    namespace standard
    {
        clock_t temp_timer;
        namespace relate
        {
            /**
             * WARNING: standard and optimized relate functions are THE SAME right now, but I keep separated because future optimization might happen.
             * 
            */

            static void APRILIntersectionFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = intersectionAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineIntersectionJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILInsideFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = insideAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineInsideJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILDisjointFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = disjointAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineDisjointJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILEqualFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = equalAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineEqualJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILMeetFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = meetAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                
                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineMeetJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILContainsFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = containsAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineContainsJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILCoversFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = coversAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                
                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineCoversJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }

            static void APRILCoveredByFilter(uint idR, uint idS) {
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    int iFilterResult = coveredByAPRILUncompressed(aprilR, aprilS);
                    if (iFilterResult != spatial_lib::INCONCLUSIVE) {
                        // true negative or true hit, return
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countAPRILResult(iFilterResult);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count result
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // refinement
                spatial_lib::time::markRefinementFilterTimer();
                spatial_lib::refineCoveredByJoin(idR, idS);
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }
        }
        
        namespace find_relation
        {
            /**
             * 
             * STANDARD APRIL (only disjoint/intersect) used for find relation queries.
             * All non-disjoint pairs are refined using DE-9IM
             * 
             */
            static void APRILFindRelation(uint idR, uint idS) {
                // printf("%u,%u\n",idR,idS);
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    // TODO: optimize (create new func that only checks AA)
                    
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                        iFilterResult = intersectionAPRILUncompressed(aprilR, aprilS);
                    } else {
                        // compressed  
                        iFilterResult = compressed::intersectionAPRIL(aprilR, aprilS);
                    }
                    
                    
                    // switch based on result
                    if (iFilterResult == spatial_lib::TRUE_NEGATIVE) {
                        // true negative, count disjoint pair
                        // measure time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        // count result
                        spatial_lib::countTopologyRelationResult(spatial_lib::TR_DISJOINT);
                        return;
                    }
                }
                // i filter ended, measure time
                spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
            
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

                // time 
                spatial_lib::time::markRefinementFilterTimer();
                // refine
                spatial_lib::refineAllRelations(idR, idS);

                // if (idR == 220594 && idS == 357429) {
                //     printf("Time passed evaluating pair %u,%u: %f seconds.\n", idR, idS, spatial_lib::time::getElapsedTime(temp_timer));
                // }


                // time
                spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            }
        }


        void IntermediateFilterEntrypoint(uint idR, uint idS) {
            // time 
            spatial_lib::time::markiFilterTimer();

            // if (idR == 220594 && idS == 357429) {
            //     temp_timer = spatial_lib::time::getNewTimer();
            // }

            // count post mbr candidate
            // spatial_lib::g_queryOutput.postMBRFilterCandidates += 1;
            // use appropriate query function
            switch (g_query->type) {
                case spatial_lib::Q_INTERSECT:
                    relate::APRILIntersectionFilter(idR, idS);
                    break;
                case spatial_lib::Q_INSIDE:
                    relate::APRILInsideFilter(idR, idS);
                    break;
                case spatial_lib::Q_DISJOINT:
                    relate::APRILDisjointFilter(idR, idS);
                    break;
                case spatial_lib::Q_EQUAL:
                    relate::APRILEqualFilter(idR, idS);
                    break;
                case spatial_lib::Q_MEET:
                    relate::APRILMeetFilter(idR, idS);
                    break;
                case spatial_lib::Q_CONTAINS:
                    relate::APRILContainsFilter(idR, idS);
                    break;
                case spatial_lib::Q_COVERS:
                    relate::APRILCoversFilter(idR, idS);
                    break;
                case spatial_lib::Q_COVERED_BY:   
                    relate::APRILCoveredByFilter(idR, idS);
                    break;
                case spatial_lib::Q_FIND_RELATION:
                    find_relation::APRILFindRelation(idR, idS);
                    break;
                default:
                    // not supported/unknown
                    printf("This combination of MBR and Intermediate filter for this query type has not been implemented yet.\n");
                    exit(-1);
                    break;
            }
        }
    }

    namespace scalability_test
    {
        namespace find_relation
        {

            namespace compressed
            {
                int intersectionAPRIL(spatial_lib::AprilDataT *aprilR, spatial_lib::AprilDataT *aprilS) {
                    // AA
                    if (!vbyte_join_compressed_sorted32(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL)) {
                        return spatial_lib::TRUE_NEGATIVE;
                    }

                    // AF
                    if (vbyte_join_compressed_sorted32(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsFULLcompressed.data(), aprilS->numIntervalsFULL)) {
                        return spatial_lib::TRUE_HIT;
                    }

                    // FA
                    if (vbyte_join_compressed_sorted32(aprilR->intervalsFULLcompressed.data(), aprilR->numIntervalsFULL, aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL)) {
                        return spatial_lib::TRUE_HIT;
                    }
                    return spatial_lib::INCONCLUSIVE;
                }

                int RinSContainmentAPRIL(spatial_lib::AprilDataT *aprilR, spatial_lib::AprilDataT *aprilS) {
                    // join AA for containment, intersection or disjoint
                    int AAresult = vbyte_join_compressed_sorted32_hybrid(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL);
                    // printf("AA: %d\n", AAresult);
                    if (AAresult == spatial_lib::IL_DISJOINT) {
                        // true negative
                        return spatial_lib::TR_DISJOINT;
                    } else if(AAresult == spatial_lib::IL_R_INSIDE_S) {
                        // all R_A intervals are inside S_A 
                        if (aprilS->numIntervalsFULL) {
                            // join AF for containment, intersection or disjoint
                            int AFresult = vbyte_join_compressed_sorted32_hybrid(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsFULLcompressed.data(), aprilS->numIntervalsFULL);
                            if (AFresult == spatial_lib::IL_R_INSIDE_S) {
                                // AF containment
                                return spatial_lib::TR_INSIDE;
                            } else if (AFresult == spatial_lib::IL_INTERSECT) {
                                // AF intersect, refinement for: inside, covered by, [true hit intersect]
                                return spatial_lib::REFINE_INSIDE_COVEREDBY_TRUEHIT_INTERSECT;
                            }
                        } else {
                            if(vbyte_join_compressed_sorted32(aprilR->intervalsFULLcompressed.data(), aprilR->numIntervalsFULL, aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL)) {
                                // intersection
                                return spatial_lib::TR_INTERSECT;
                            }
                        }
                    } else {
                        // intersection
                        if(vbyte_join_compressed_sorted32(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsFULLcompressed.data(), aprilS->numIntervalsFULL)){
                            // intersection
                            return spatial_lib::TR_INTERSECT;
                        } else if(vbyte_join_compressed_sorted32(aprilR->intervalsFULLcompressed.data(), aprilR->numIntervalsFULL, aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL)) {
                            // intersection
                            return spatial_lib::TR_INTERSECT;
                        }
                    }

                    // refinement for: disjoint, inside, covered by, meet, intersect
                    return spatial_lib::REFINE_DISJOINT_INSIDE_COVEREDBY_MEET_INTERSECT;
                }

                int SinRContainmentAPRIL(spatial_lib::AprilDataT *aprilR, spatial_lib::AprilDataT *aprilS) {
                    // join AA for containment, intersection or disjoint
                    int AAresult = vbyte_join_compressed_sorted32_hybrid(aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL, aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL);
                    if (AAresult == spatial_lib::IL_DISJOINT) {
                        // true negative
                        return spatial_lib::TR_DISJOINT;
                    } else if(AAresult == spatial_lib::IL_R_INSIDE_S) {
                        // all S_A intervals are inside R_A 
                        if (aprilR->numIntervalsFULL) {
                            // join AF for containment, intersection or disjoint
                            int AFresult = vbyte_join_compressed_sorted32_hybrid(aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL, aprilR->intervalsFULLcompressed.data(), aprilR->numIntervalsFULL);
                            if (AFresult == spatial_lib::IL_R_INSIDE_S) {
                                // AF containment (dont confuse IL_R_INSIDE_S, it is based on how you pass arguments to the hybrid function)
                                return spatial_lib::TR_CONTAINS;
                            } else if (AFresult == spatial_lib::IL_INTERSECT) {
                                // AF intersect, refinement for: contains, covers, [true hit intersect]
                                return spatial_lib::REFINE_CONTAINS_COVERS_TRUEHIT_INTERSECT;
                            }
                        } else {
                            if(vbyte_join_compressed_sorted32(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsFULLcompressed.data(), aprilS->numIntervalsFULL)) {
                                // check AF intersection
                                return spatial_lib::TR_INTERSECT;
                            }
                        }
                    } else {
                        // intersection
                        if(vbyte_join_compressed_sorted32(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsFULLcompressed.data(), aprilS->numIntervalsFULL)){
                            // intersection
                            return spatial_lib::TR_INTERSECT;
                        } else if(vbyte_join_compressed_sorted32(aprilR->intervalsFULLcompressed.data(), aprilR->numIntervalsFULL, aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL)) {
                            // intersection
                            return spatial_lib::TR_INTERSECT;
                        }
                    }

                    // refinement for: disjoint, inside, covered by, meet, intersect
                    return spatial_lib::REFINE_DISJOINT_CONTAINS_COVERS_MEET_INTERSECT;
                }

                int MBRIntersectionAPRIL(spatial_lib::AprilDataT *aprilR, spatial_lib::AprilDataT *aprilS) {
                    // use regular APRIL
                    int aprilResult = intersectionAPRIL(aprilR, aprilS);
                    switch(aprilResult) {
                        case spatial_lib::TRUE_NEGATIVE:
                            return spatial_lib::TR_DISJOINT;
                        case spatial_lib::TRUE_HIT:
                            return spatial_lib::TR_INTERSECT;
                        case spatial_lib::INCONCLUSIVE:
                            return spatial_lib::REFINE_DISJOINT_MEET_INTERSECT;
                    }
                }

                int equalMBRsAPRIL(uint idR, uint idS, spatial_lib::AprilDataT *aprilR, spatial_lib::AprilDataT *aprilS) {
                    // printf("%u,%u\n", idR, idS);
                    // AA join to look for exact relationship between the lists
                    // int AAresult = joinIntervalListsSymmetricalOptimized(aprilR->intervalsALL, aprilR->numIntervalsALL, aprilS->intervalsALL, aprilS->numIntervalsALL);
                    int AAresult = vbyte_join_compressed_sorted32_symmetrical(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL);
                    if (AAresult == spatial_lib::IL_MATCH) {
                        // refine for equal, covered by, covers and true hit intersect
                        return spatial_lib::REFINE_EQUAL_COVERS_COVEREDBY_TRUEHIT_INTERSECT;
                    } else if (AAresult == spatial_lib::IL_R_INSIDE_S) {
                        int AFresult = vbyte_join_compressed_sorted32_hybrid(aprilR->intervalsALLcompressed.data(), aprilR->numIntervalsALL, aprilS->intervalsFULLcompressed.data(), aprilS->numIntervalsFULL);
                        // printf("AF: %d\n", AFresult);
                        if (AFresult == spatial_lib::IL_R_INSIDE_S) {
                            // true hit covered by (return inside because reasons)
                            // return spatial_lib::TR_COVERED_BY;
                            return spatial_lib::TR_INSIDE;
                        } 
                        // AA no containment, AF disjoint or intersect
                        return spatial_lib::REFINE_COVEREDBY_TRUEHIT_INTERSECT;
                    } else if(AAresult == spatial_lib::IL_S_INSIDE_R) {
                        int FAresult = vbyte_join_compressed_sorted32_hybrid(aprilS->intervalsALLcompressed.data(), aprilS->numIntervalsALL, aprilR->intervalsFULLcompressed.data(), aprilR->numIntervalsFULL);
                        // printf("FA: %d\n", FAresult);
                        // in this case R is S because joinIntervalsHybrid handles the first list as R and the second as S
                        // and only checks for assymetrical containment R in S
                        if (FAresult == spatial_lib::IL_R_INSIDE_S) {
                            // FA true hit covers, return contains
                            // return spatial_lib::TR_COVERS;
                            return spatial_lib::TR_CONTAINS;
                        } 
                        // AA no containment, FA disjoint or intersect
                        return spatial_lib::REFINE_COVERS_TRUEHIT_INTERSECT;
                    } else {
                        // AA no containment, true hit intersect because equal MBRs
                        // maybe a special case of meet, so refine first
                        if(spatial_lib::isMeet(idR, idS)){
                            return spatial_lib::TR_MEET;
                        }
                        return spatial_lib::TR_INTERSECT;
                    }
                }
            }

            /************************************************/
            /**
             * SPECIALIZED TOPOLOGY FOR SCALABILITY TESTS
            */
            /************************************************/
            static void specializedTopologyRinSContainment(uint idR, uint idS) {
                // get bucket ID
                uint bucketID = spatial_lib::getBucketOfPair(idR, idS);

                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                        iFilterResult = RinSContainmentAPRILUncompressed(aprilR, aprilS);
                    } else {
                        // compressed
                        iFilterResult = compressed::RinSContainmentAPRIL(aprilR, aprilS);
                    }
                    
                    
                    
                    
                    // switch based on result
                    switch(iFilterResult) {
                        // true hits, count and return
                        case spatial_lib::TR_INSIDE:
                        case spatial_lib::TR_DISJOINT:
                        case spatial_lib::TR_INTERSECT:
                            // result
                            spatial_lib::countTopologyRelationResult(iFilterResult);

                            // if(iFilterResult == spatial_lib::TR_INTERSECT){
                            //     printf("%u,%u\n", idR, idS);
                            // }

                            // time
                            double elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                            spatial_lib::g_queryOutput.iFilterTime += elapsedTime;
                            spatial_lib::addIFilterTimeToBucket(bucketID, elapsedTime);
                            return;
                    }
                }
                // i filter ended, measure time
                double elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                spatial_lib::g_queryOutput.iFilterTime += elapsedTime;
                spatial_lib::addIFilterTimeToBucket(bucketID, elapsedTime);
            
                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                spatial_lib::countInconclusiveToBucket(bucketID);

                int relation;
                // switch based on result
                switch(iFilterResult) {            
                    // inconclusive, do selective refinement
                    // result holds what type of refinement needs to be made
                    case spatial_lib::REFINE_INSIDE_COVEREDBY_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineInsideCoveredbyTruehitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_DISJOINT_INSIDE_COVEREDBY_MEET_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineDisjointInsideCoveredbyMeetIntersect(idR, idS);
                        break;
                    
                }
                // time
                elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
                spatial_lib::g_queryOutput.refinementTime += elapsedTime;
                spatial_lib::addRefinementTimeToBucket(bucketID, elapsedTime);
                
                // count the result
                spatial_lib::countTopologyRelationResult(relation);
                // if(relation == spatial_lib::TR_INTERSECT){
                //     printf("%u,%u\n", idR, idS);
                // }
            }

            static void specializedTopologySinRContainment(uint idR, uint idS) {
                // get bucket ID
                uint bucketID = spatial_lib::getBucketOfPair(idR, idS);
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                        iFilterResult = SinRContainmentAPRILUncompressed(aprilR, aprilS);
                    } else {
                        // compressed
                        iFilterResult = compressed::SinRContainmentAPRIL(aprilR, aprilS);
                    }

                    // switch based on result
                    switch(iFilterResult) {
                        // true hits, count and return
                        case spatial_lib::TR_CONTAINS:
                        case spatial_lib::TR_DISJOINT:
                        case spatial_lib::TR_INTERSECT:
                            // if(iFilterResult == spatial_lib::TR_INTERSECT){
                            //     printf("%u,%u\n", idR, idS);
                            // }
                            // result
                            spatial_lib::countTopologyRelationResult(iFilterResult);
                            // time
                            double elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                            spatial_lib::g_queryOutput.iFilterTime += elapsedTime;
                            spatial_lib::addIFilterTimeToBucket(bucketID, elapsedTime);
                            return;
                    }
                }
                // i filter ended, measure time
                double elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                spatial_lib::g_queryOutput.iFilterTime += elapsedTime;
                spatial_lib::addIFilterTimeToBucket(bucketID, elapsedTime);

                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                spatial_lib::countInconclusiveToBucket(bucketID);

                int relation;
                // switch based on result
                switch(iFilterResult) {            
                    // inconclusive, do selective refinement
                    // result holds what type of refinement needs to be made
                    case spatial_lib::REFINE_CONTAINS_COVERS_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineContainsCoversTruehitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_DISJOINT_CONTAINS_COVERS_MEET_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineDisjointContainsCoversMeetIntersect(idR, idS);
                        break;
                    
                }
                // time
                elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
                spatial_lib::g_queryOutput.refinementTime += elapsedTime;
                spatial_lib::addRefinementTimeToBucket(bucketID, elapsedTime);
                
                // count the result
                spatial_lib::countTopologyRelationResult(relation);
                // if(relation == spatial_lib::TR_INTERSECT){
                //     printf("%u,%u\n", idR, idS);
                // }
            }

            static void specializedTopologyIntersection(uint idR, uint idS) {
                // get bucket ID
                uint bucketID = spatial_lib::getBucketOfPair(idR, idS);

                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                        iFilterResult = MBRIntersectionAPRILUncompressed(aprilR, aprilS);
                    } else {
                        // compressed
                        iFilterResult = compressed::MBRIntersectionAPRIL(aprilR, aprilS);
                    }

                    // switch based on result
                    switch(iFilterResult) {
                        // true hits, count and return
                        case spatial_lib::TR_DISJOINT:
                        case spatial_lib::TR_INTERSECT:
                            // if(iFilterResult == spatial_lib::TR_INTERSECT){
                            //     printf("%u,%u\n", idR, idS);
                            // }
                            // result
                            spatial_lib::countTopologyRelationResult(iFilterResult);
                            // time
                            double elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                            spatial_lib::g_queryOutput.iFilterTime += elapsedTime;
                            spatial_lib::addIFilterTimeToBucket(bucketID, elapsedTime);
                            return;
                    }
                }
                // i filter ended, measure time
                double elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                spatial_lib::g_queryOutput.iFilterTime += elapsedTime;
                spatial_lib::addIFilterTimeToBucket(bucketID, elapsedTime);

                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                spatial_lib::countInconclusiveToBucket(bucketID);

                // time 
                spatial_lib::time::markRefinementFilterTimer();
                // refine
                int relation = spatial_lib::refineDisjointMeetIntersect(idR, idS);
                // time
                elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
                spatial_lib::g_queryOutput.refinementTime += elapsedTime;
                spatial_lib::addRefinementTimeToBucket(bucketID, elapsedTime);
                // count the result
                spatial_lib::countTopologyRelationResult(relation);
                // if(relation == spatial_lib::TR_INTERSECT){
                //     printf("%u,%u\n", idR, idS);
                // }
            }

            static void specializedTopologyEqual(uint idR, uint idS) {
                // get bucket ID
                uint bucketID = spatial_lib::getBucketOfPair(idR, idS);
                // get common sections
                std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
                // for each common section
                int iFilterResult;
                for (auto &sectionID : commonSections) {
                    // fetch the APRIL of R and S for this section
                    spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                    spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);

                    // use APRIL intermediate filter
                    if (g_query->R.aprilConfig.compression == spatial_lib::C_COMPRESSION_DISABLED) {
                        // uncompressed
                    iFilterResult = equalMBRsAPRILUncompressed(idR, idS, aprilR, aprilS);
                    } else {
                        // compressed
                        iFilterResult = compressed::equalMBRsAPRIL(idR, idS, aprilR, aprilS);
                    }
                
                    // switch based on result
                    switch(iFilterResult) {
                        // true hits, count and return
                        case spatial_lib::TR_DISJOINT:
                        case spatial_lib::TR_EQUAL:
                        case spatial_lib::TR_COVERED_BY:
                        case spatial_lib::TR_COVERS:
                        case spatial_lib::TR_INTERSECT:
                            // if(iFilterResult == spatial_lib::TR_INTERSECT) {
                            //     printf("%u,%u\n",idR,idS);
                            // }
                            // result
                            spatial_lib::countTopologyRelationResult(iFilterResult);
                            // time
                            double elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                            spatial_lib::g_queryOutput.iFilterTime += elapsedTime;
                            spatial_lib::addIFilterTimeToBucket(bucketID, elapsedTime);
                            return;
                    }
                }
                // i filter ended, measure time
                double elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                spatial_lib::g_queryOutput.iFilterTime += elapsedTime;
                spatial_lib::addIFilterTimeToBucket(bucketID, elapsedTime);

                // count refinement candidate
                spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);
                spatial_lib::countInconclusiveToBucket(bucketID);

                int relation;
                // switch based on result
                switch(iFilterResult) {            
                    // inconclusive, do selective refinement
                    // result holds what type of refinement needs to be made
                    case spatial_lib::REFINE_COVEREDBY_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineCoveredbyTrueHitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_COVERS_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineCoversTrueHitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_COVERS_COVEREDBY_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineCoversCoveredByTrueHitIntersect(idR, idS);
                        break;
                    case spatial_lib::REFINE_EQUAL_COVERS_COVEREDBY_TRUEHIT_INTERSECT:
                        // time 
                        spatial_lib::time::markRefinementFilterTimer();
                        // refine
                        relation = spatial_lib::refineEqualCoversCoveredByTrueHitIntersect(idR, idS);
                        break;
                    
                }
                // time
                elapsedTime = spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
                spatial_lib::g_queryOutput.refinementTime += elapsedTime;
                spatial_lib::addRefinementTimeToBucket(bucketID, elapsedTime);
                // count the result
                spatial_lib::countTopologyRelationResult(relation);
                // if(relation == spatial_lib::TR_INTERSECT) {
                //     printf("%u,%u\n",idR,idS);
                // }
            }
        }

        void IntermediateFilterEntrypoint(uint idR, uint idS, int relationCase) {
            // time 
            spatial_lib::time::markiFilterTimer();

            // count post mbr candidate
            // spatial_lib::g_queryOutput.postMBRFilterCandidates += 1;
            // use appropriate query function
            switch (relationCase) {
                case spatial_lib::MBR_R_IN_S:
                    find_relation::specializedTopologyRinSContainment(idR, idS);
                    break;
                case spatial_lib::MBR_S_IN_R:
                    find_relation::specializedTopologySinRContainment(idR, idS);
                    break;
                case spatial_lib::MBR_EQUAL:
                    find_relation::specializedTopologyEqual(idR, idS);
                    break;
                case spatial_lib::MBR_INTERSECT:
                    find_relation::specializedTopologyIntersection(idR, idS);
                    break;
                default:
                    // not supported/unknown
                    exit(-1);
                    break;
            }
        }   
    }

    namespace on_the_fly
    {

        /************************************************/
        /**
         * ON THE FLY RASTERIZATION WITH SPECIALIZED 
         * TOPOLOGY FROM MBR FILTER
        */
        /************************************************/
        static spatial_lib::AprilDataT* createAPRILForObject(uint recID, uint sectionID, spatial_lib::DatasetT *dataset, bool R) {
            // load polygon
            spatial_lib::bg_polygon polygon = spatial_lib::loadBoostPolygonByIDandFlag(recID, R);

            // get section by ID
            spatial_lib::SectionT* sec = spatial_lib::getSectionByID(*dataset, sectionID);

            // set rasterize dataspace
            rasterizerlib::setDataspace(sec->rasterxMin, sec->rasteryMin, sec->rasterxMax, sec->rasteryMax);

            // rasterize polygon
            spatial_lib::AprilDataT aprilData = rasterizerlib::generateAPRILForBoostGeometry(polygon);

            // store to section 
            spatial_lib::addAprilDataToApproximationDataMap(*dataset, sectionID, recID, &aprilData);
            // stored, now delete allocation done by the rasterizer

            // count rasterization
            spatial_lib::g_queryOutput.rasterizationsDone += 1;

            // return object (from map, local reference obsolete after store)
            return spatial_lib::getAprilDataBySectionAndObjectIDs(*dataset, sectionID, recID);
        }

        static void specializedTopologyRinSContainmentOTF(uint idR, uint idS) {
            // get common sections
            std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
            // for each common section
            int iFilterResult;
            for (auto &sectionID : commonSections) {
                // check if APRIL has been created for these objects in this section
                spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);
                if (aprilR == nullptr) {
                    // create and store APRIL for idR
                    aprilR = createAPRILForObject(idR, sectionID, &g_query->R, true);
                }
                if (aprilS == nullptr) {
                    // create and store APRIL for idS
                    aprilS = createAPRILForObject(idS, sectionID, &g_query->S, false);
                }

                // use FIND RELATION APRIL intermediate filter normally
                iFilterResult = RinSContainmentAPRILUncompressed(aprilR, aprilS);
                // switch based on result
                switch(iFilterResult) {
                    // true hits, count and return
                    case spatial_lib::TR_INSIDE:
                    case spatial_lib::TR_DISJOINT:
                    case spatial_lib::TR_INTERSECT:
                        // result
                        spatial_lib::countTopologyRelationResult(iFilterResult);

                        // if(iFilterResult == spatial_lib::TR_INTERSECT){
                        //     printf("%u,%u\n", idR, idS);
                        // }

                        // time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        return;
                }
            }
            // i filter ended, measure time
            spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
        
            // count refinement candidate
            spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

            int relation;
            // switch based on result
            switch(iFilterResult) {            
                // inconclusive, do selective refinement
                // result holds what type of refinement needs to be made
                case spatial_lib::REFINE_INSIDE_COVEREDBY_TRUEHIT_INTERSECT:
                    // time 
                    spatial_lib::time::markRefinementFilterTimer();
                    // refine
                    relation = spatial_lib::refineInsideCoveredbyTruehitIntersect(idR, idS);
                    break;
                case spatial_lib::REFINE_DISJOINT_INSIDE_COVEREDBY_MEET_INTERSECT:
                    // time 
                    spatial_lib::time::markRefinementFilterTimer();
                    // refine
                    relation = spatial_lib::refineDisjointInsideCoveredbyMeetIntersect(idR, idS);
                    break;
                
            }
            // time
            spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            // count the result
            spatial_lib::countTopologyRelationResult(relation);
            // if(relation == spatial_lib::TR_INTERSECT){
            //     printf("%u,%u\n", idR, idS);
            // }
        }

        static void specializedTopologySinRContainmentOTF(uint idR, uint idS) {
            // get common sections
            std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
            // for each common section
            int iFilterResult;
            for (auto &sectionID : commonSections) {
                // check if APRIL has been created for these objects in this section
                spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);
                
                if (aprilR == nullptr) {
                    // create and store APRIL for idR
                    aprilR = createAPRILForObject(idR, sectionID, &g_query->R, true);
                }
                if (aprilS == nullptr) {
                    // create and store APRIL for idS
                    aprilS = createAPRILForObject(idS, sectionID, &g_query->S, false);
                }

                // use APRIL intermediate filter
                iFilterResult = SinRContainmentAPRILUncompressed(aprilR, aprilS);
                // switch based on result
                switch(iFilterResult) {
                    // true hits, count and return
                    case spatial_lib::TR_CONTAINS:
                    case spatial_lib::TR_DISJOINT:
                    case spatial_lib::TR_INTERSECT:
                        // if(iFilterResult == spatial_lib::TR_INTERSECT){
                        //     printf("%u,%u\n", idR, idS);
                        // }
                        // result
                        spatial_lib::countTopologyRelationResult(iFilterResult);
                        // time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        return;
                }
            }
            // i filter ended, measure time
            spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
        
            // count refinement candidate
            spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

            int relation;
            // switch based on result
            switch(iFilterResult) {            
                // inconclusive, do selective refinement
                // result holds what type of refinement needs to be made
                case spatial_lib::REFINE_CONTAINS_COVERS_TRUEHIT_INTERSECT:
                    // time 
                    spatial_lib::time::markRefinementFilterTimer();
                    // refine
                    relation = spatial_lib::refineContainsCoversTruehitIntersect(idR, idS);
                    break;
                case spatial_lib::REFINE_DISJOINT_CONTAINS_COVERS_MEET_INTERSECT:
                    // time 
                    spatial_lib::time::markRefinementFilterTimer();
                    // refine
                    relation = spatial_lib::refineDisjointContainsCoversMeetIntersect(idR, idS);
                    break;
                
            }
            // time
            spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            // count the result
            spatial_lib::countTopologyRelationResult(relation);
            // if(relation == spatial_lib::TR_INTERSECT){
            //     printf("%u,%u\n", idR, idS);
            // }
        }

        static void specializedTopologyEqualOTF(uint idR, uint idS) {
            // get common sections
            std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
            // for each common section
            int iFilterResult;
            for (auto &sectionID : commonSections) {
                // check if APRIL has been created for these objects in this section
                spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);
                if (aprilR == nullptr) {
                    // create and store APRIL for idR
                    aprilR = createAPRILForObject(idR, sectionID, &g_query->R, true);
                }
                if (aprilS == nullptr) {
                    // create and store APRIL for idS
                    aprilS = createAPRILForObject(idS, sectionID, &g_query->S, false);
                }

                // use APRIL intermediate filter
                iFilterResult = equalMBRsAPRILUncompressed(idR, idS, aprilR, aprilS);
            
                // switch based on result
                switch(iFilterResult) {
                    // true hits, count and return
                    case spatial_lib::TR_DISJOINT:
                    case spatial_lib::TR_EQUAL:
                    case spatial_lib::TR_COVERED_BY:
                    case spatial_lib::TR_COVERS:
                    case spatial_lib::TR_INTERSECT:
                        // if(iFilterResult == spatial_lib::TR_INTERSECT) {
                        //     printf("%u,%u\n",idR,idS);
                        // }
                        // result
                        spatial_lib::countTopologyRelationResult(iFilterResult);
                        // time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        return;
                }
            }
            // i filter ended, measure time
            spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
        
            // count refinement candidate
            spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

            int relation;
            // switch based on result
            switch(iFilterResult) {            
                // inconclusive, do selective refinement
                // result holds what type of refinement needs to be made
                case spatial_lib::REFINE_COVEREDBY_TRUEHIT_INTERSECT:
                    // time 
                    spatial_lib::time::markRefinementFilterTimer();
                    // refine
                    relation = spatial_lib::refineCoveredbyTrueHitIntersect(idR, idS);
                    break;
                case spatial_lib::REFINE_COVERS_TRUEHIT_INTERSECT:
                    // time 
                    spatial_lib::time::markRefinementFilterTimer();
                    // refine
                    relation = spatial_lib::refineCoversTrueHitIntersect(idR, idS);
                    break;
                case spatial_lib::REFINE_COVERS_COVEREDBY_TRUEHIT_INTERSECT:
                    // time 
                    spatial_lib::time::markRefinementFilterTimer();
                    // refine
                    relation = spatial_lib::refineCoversCoveredByTrueHitIntersect(idR, idS);
                    break;
                case spatial_lib::REFINE_EQUAL_COVERS_COVEREDBY_TRUEHIT_INTERSECT:
                    // time 
                    spatial_lib::time::markRefinementFilterTimer();
                    // refine
                    relation = spatial_lib::refineEqualCoversCoveredByTrueHitIntersect(idR, idS);
                    break;
            }
            // time
            spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            // count the result
            spatial_lib::countTopologyRelationResult(relation);
            // if(relation == spatial_lib::TR_INTERSECT) {
            //     printf("%u,%u\n",idR,idS);
            // }
        }

        static void specializedTopologyIntersectionOTF(uint idR, uint idS) {
            // get common sections
            std::vector<uint> commonSections = spatial_lib::getCommonSectionIDsOfObjects(g_query->R, g_query->S, idR, idS);
            // for each common section
            int iFilterResult;
            for (auto &sectionID : commonSections) {
                // check if APRIL has been created for these objects in this section
                spatial_lib::AprilDataT* aprilR = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->R, sectionID, idR);
                spatial_lib::AprilDataT* aprilS = spatial_lib::getAprilDataBySectionAndObjectIDs(g_query->S, sectionID, idS);
                if (aprilR == nullptr) {
                    // create and store APRIL for idR
                    aprilR = createAPRILForObject(idR, sectionID, &g_query->R, true);
                }
                if (aprilS == nullptr) {
                    // create and store APRIL for idS
                    aprilS = createAPRILForObject(idS, sectionID, &g_query->S, false);
                }

                // use APRIL intermediate filter
                iFilterResult = MBRIntersectionAPRILUncompressed(aprilR, aprilS);
                // switch based on result
                switch(iFilterResult) {
                    // true hits, count and return
                    case spatial_lib::TR_DISJOINT:
                    case spatial_lib::TR_INTERSECT:
                        // if(iFilterResult == spatial_lib::TR_INTERSECT){
                        //     printf("%u,%u\n", idR, idS);
                        // }
                        // result
                        spatial_lib::countTopologyRelationResult(iFilterResult);
                        // time
                        spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
                        return;
                }
            }
            // i filter ended, measure time
            spatial_lib::g_queryOutput.iFilterTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.iFilterTimer);
        
            // count refinement candidate
            spatial_lib::countAPRILResult(spatial_lib::INCONCLUSIVE);

            // time 
            spatial_lib::time::markRefinementFilterTimer();
            // refine
            int relation = spatial_lib::refineDisjointMeetIntersect(idR, idS);
            // time
            spatial_lib::g_queryOutput.refinementTime += spatial_lib::time::getElapsedTime(spatial_lib::time::g_timer.refTimer);
            // count the result
            spatial_lib::countTopologyRelationResult(relation);
            // if(relation == spatial_lib::TR_INTERSECT){
            //     printf("%u,%u\n", idR, idS);
            // }
        }

        void IntermediateFilterEntrypointOTF(uint idR, uint idS, int relationCase) {
            // time 
            spatial_lib::time::markiFilterTimer();

            // count post mbr candidate
            // spatial_lib::g_queryOutput.postMBRFilterCandidates += 1;
            // use appropriate query function
            switch (relationCase) {
                case spatial_lib::MBR_R_IN_S:
                    specializedTopologyRinSContainmentOTF(idR, idS);
                    break;
                case spatial_lib::MBR_S_IN_R:
                    specializedTopologySinRContainmentOTF(idR, idS);
                    break;
                case spatial_lib::MBR_EQUAL:
                    specializedTopologyEqualOTF(idR, idS);
                    break;
                case spatial_lib::MBR_INTERSECT:
                    specializedTopologyIntersectionOTF(idR, idS);
                    break;
                default:
                    // not supported/unknown
                    exit(-1);
                    break;
            }
            

        }
    }

}