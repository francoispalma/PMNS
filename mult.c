#include <stdio.h>
#include <stdint.h>

#include "montgom.h"


int main(void)
{
	poly A, B, C;
	init_polys(5, &A, &B, NULL);
	init_poly(10, &C);
	
	const char a[] = "163dbed3b4fbeee3bb542bc62983a51f4ecb077c3af9a6d451e1c4b6cd0a99563d55",
	b[] = "c64be1b07fc889170737a3bbd0501940eb5cdffbb228d09f6c4527812aa64b8cde15";
	
	convert_string_to_poly(&A, a);
	print(A);
	convert_string_to_poly(&B, b);
	print(B);
	
	free_polys(A, B, C, NULL);
	return 0;
}
