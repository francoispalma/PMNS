#ifndef PMNS_H
#define PMNS_H

void randpoly(poly P);
void pmns_montg_mult(poly res, const poly A,
	const poly B);
void pmns_rtl_sqandmult(poly res, const poly base,
	const mpnum exponent);
void pmns_ltr_sqandmult(poly res, const poly base,
	const mpnum exponent);
void pmns_montg_ladder(poly res, const poly base,
	const mpnum exponent);
void pmns_sqandmult(poly res, const char* base,
	const char* exponent);
void convert_binary_to_pmns(poly res, const mpnum op);
void convert_string_to_pmns(poly res, const char* string);
void convert_pmns_to_binary(mpnum* res, const poly P);

#endif
