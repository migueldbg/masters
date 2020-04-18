#include <TCut.h>
#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

/* *********************************************************************************************************************** *
 * File: RunsCheck.C
 *
 * Author: Miguel Del Ben Galdiano.
 * Date of Creation : April 18 2020.
 *
 * Summary of File:
 *
 *    The goal of this macro is to run a preliminary analysis to a set of runs (from 1489 up to 1521). It will also serve
 *    as a learning experience to figure out how to join TTrees with the same structure that come from different root
 *    files. The runs will be separated in two groups: those with drift field turned off (1489 to 1500) and those with
 *    drift field turned on (1501 to 1521). The first group will be further divided based on the DAQ configuration, while
 *    the second will be divide based on the CH2 target density. I believe these groups will show similar behaviour and
 *    the subdivisions will prove to be irrelevant, but nonetheless.
 *
 *    TABLE OF DIVISIONS:
 *
 *    ::                      FIELD ON                       ::                   FIELD OFF                    ::
 *    ::  SiMaster_TPCslave   Si_Standalone    Si_and_TPC    ::    401 ug/cm2     426 ug/cm2     423 ug/cm2    ::
 *    ::    1489 to 1494      1495 to 1496    1499 to 1500   ::   1501 to 1510   1511 to 1516   1519 to 1521   ::
 *
 *    DISCARDED RUNS:
 *      > 1497: Au (90 um/cm2) run.
 *      > 1498: Faulty channels, told not to use.
 *      > 1505: Au (90 um/cm2) run.
 *      > 1507: JUNK.
 *      > 1517: JUNK.
 *      > 1519: JUNK.
 *
 *    SPECIFIC RUN NOTES:
 *      > 1501: Rise field! Beam not stable.
 *
 * *********************************************************************************************************************** */

 
