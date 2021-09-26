#include "gf_cal.h"
#include "ee_decoding.h"
#include "encoding.h"

void main()
{
	unsigned char i = 0;

	systematic_encoding();

	/*transmission through channel*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		received_polynomial[i] = gf_add(encoded_polynomial[i], error_polynomial[i]);
	}
	
	ee_decoding();

	return;
}
