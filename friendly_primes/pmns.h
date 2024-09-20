#ifndef PMNS_H_INCLUDED
#define PMNS_H_INCLUDED

#include <stdint.h>
#include "structs.h"

extern uint8_t ALPHA;

void setpmnsparams(uint16_t newN, uint64_t newGAMMA, int16_t newLAMBDA);
void setpmnsalpha(uint8_t newALPHA);
void setpmnspsi(uint64_t psi);
void setpmnsaux(uint64_t newGIMEL, int16_t newLAMED);
void setpmnslastcol(uint64_t newcol[]);
void setpmnssparseinv(uint64_t newGL1, uint64_t newGL2, uint64_t newOL1);

void pmns_linearred_mult(restrict poly res, const restrict poly A, const restrict poly B);
void pmns_doublesparse_mult(restrict poly res, const restrict poly A, const restrict poly B);

void pmns_linearred_mult9(restrict poly res, const restrict poly A, const restrict poly B);
void pmns_doublesparse_mult9(restrict poly res, const restrict poly A, const restrict poly B);

void pmns_linearred_multa9(restrict poly res, const restrict poly A, const restrict poly B);
void pmns_doublesparse_multa9(restrict poly res, const restrict poly A, const restrict poly B);

#endif
