#include <stdio.h>
#include <string.h>
#include "gf_cal.h"
#include "ee_decoding.h"

unsigned char received_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

unsigned char output_polynomial[CODEWORD_LEN] =
{
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF,
	0xFF
};

int ee_decoding()
{
	unsigned char i = 0, j = 0, tmp = 0xFF, tmp_sum = 0xFF, lambda_root = 0;
	unsigned char syndrome[CODEWORD_LEN - MESSAGE_LEN];
	unsigned char lambda[(CODEWORD_LEN - MESSAGE_LEN) / 2 + 1];
	unsigned char omega[(CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN)];
	unsigned char err_location[CODEWORD_LEN];
	unsigned char lambda_dev[(CODEWORD_LEN - MESSAGE_LEN) / 2];
	unsigned char err_mag[CODEWORD_LEN];
	unsigned char codeword[CODEWORD_LEN];
	unsigned char r_y[CODEWORD_LEN - MESSAGE_LEN + 1];
	unsigned char r_z[CODEWORD_LEN - MESSAGE_LEN + 1];
	unsigned char r_a[CODEWORD_LEN - MESSAGE_LEN + 1];
	unsigned char q_a[CODEWORD_LEN - MESSAGE_LEN + 1];
	#if 0
	unsigned char s_y[CODEWORD_LEN - MESSAGE_LEN + 1];
	unsigned char s_z[CODEWORD_LEN - MESSAGE_LEN + 1];
	unsigned char s_a[CODEWORD_LEN - MESSAGE_LEN + 1];
#endif	
	unsigned char t_y[CODEWORD_LEN - MESSAGE_LEN + 1];
	unsigned char t_z[CODEWORD_LEN - MESSAGE_LEN + 1];
	unsigned char t_a[CODEWORD_LEN - MESSAGE_LEN + 1];
	unsigned char product_tmp[CODEWORD_LEN - MESSAGE_LEN + 1];
	
	memset(syndrome, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));
	memset(lambda, 0x0, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + 1));
	memset(omega, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN)));
	memset(err_location, 0x0, sizeof(unsigned char) * CODEWORD_LEN);
	memset(lambda_dev, 0xFF, sizeof(unsigned char) * ((CODEWORD_LEN - MESSAGE_LEN) / 2));
	memset(err_mag, 0xFF, sizeof(unsigned char) * CODEWORD_LEN);
	memcpy(codeword, received_polynomial, sizeof(unsigned char) * CODEWORD_LEN);

	printf("Received:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", received_polynomial[i]);
	}
	printf("\n");

	/*compute the syndrome polynomial from received symbols*/
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		tmp = 0xFF;
		tmp_sum = 0xFF;
		for(j = 0; j < CODEWORD_LEN; j++)
		{
			tmp = gf_multp(received_polynomial[j], (i + 1) * j);
			tmp_sum = gf_add(tmp, tmp_sum);
			//printf("%x %x\n", tmp, tmp_sum);
		}
		syndrome[i] = tmp_sum;
	}
	printf("Syndrome Polynomial:\n");
	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		printf("%x ", syndrome[i]);
	}
	printf("\n");

	memset(r_a, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
	//memset(s_a, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
	memset(t_a, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
	
	memset(q_a, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
	memset(product_tmp, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));

	memcpy(r_z, syndrome, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));//no +1
	r_z[CODEWORD_LEN - MESSAGE_LEN] = 0xFF;
	memset(r_y, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN));//no +1
	r_y[CODEWORD_LEN - MESSAGE_LEN] = 0;

#if 0
	memset(s_z, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));//no +1
	memset(s_y, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));//no +1
	s_y[0] = 0;
#endif

	memset(t_z, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));//no +1
	t_z[0] = 0;
	memset(t_y, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));//no +1

	for(i = 0; i < CODEWORD_LEN - MESSAGE_LEN; i++)
	{
		gf_div_q_r(r_y, CODEWORD_LEN - MESSAGE_LEN + 1,
					r_z, CODEWORD_LEN - MESSAGE_LEN + 1,
					q_a, CODEWORD_LEN - MESSAGE_LEN + 1,
					r_a, CODEWORD_LEN - MESSAGE_LEN + 1);

		gf_multp_poly(t_z, CODEWORD_LEN - MESSAGE_LEN + 1,
					   q_a, CODEWORD_LEN - MESSAGE_LEN + 1,
					   product_tmp, CODEWORD_LEN - MESSAGE_LEN + 1);
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN + 1; j++)
		{
			t_a[j] = gf_add(t_y[j], product_tmp[j]);
			//printf("%x %x %x\n", t_y[j], product_tmp[j], t_a[j]);
		}

		printf("Q_A:\n");
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN + 1; j++)
		{
			printf("%x ", q_a[j]);
		}
		printf("\n");
		printf("R_A:\n");
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN + 1; j++)
		{
			printf("%x ", r_a[j]);
		}
		printf("\n");
		printf("T_A:\n");
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN + 1; j++)
		{
			printf("%x ", t_a[j]);
		}
		printf("\n");

#if 0
		tmp = 0;
		for(j = 0; j < CODEWORD_LEN - MESSAGE_LEN + 1; j++)
		{
			if(0xFF != r_a[j])
			{
				tmp = tmp + 1;
			}
		}
		if(0 == tmp)
		{
			printf("decoded OK: %d\n", i);
		}
#endif
		if(((CODEWORD_LEN - MESSAGE_LEN) / 2) > gf_degree(r_a, CODEWORD_LEN - MESSAGE_LEN + 1))
		{
			printf("Iteration OK: %d\n", i);
			break;
		}

		memset(q_a, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
		memset(product_tmp, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
		
		memcpy(r_y, r_z, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
		memcpy(r_z, r_a, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
		memset(r_a, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
#if 0
		memcpy(s_y, s_z, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
		memcpy(s_z, s_a, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
		memset(s_a, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
#endif
		memcpy(t_y, t_z, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
		memcpy(t_z, t_a, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
		memset(t_a, 0xFF, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN + 1));
	}

	memcpy(omega, r_a, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN));
	memcpy(lambda, t_a, sizeof(unsigned char) * (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1);

	printf("Lambda Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; i++)
	{
		printf("%x ", lambda[i]);
	}
	printf("\n");
	printf("Omega Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN); i++)
	{
		printf("%x ", omega[i]);
	}
	printf("\n");

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		lambda_root = 0xFF;
		for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + 1; j++)
		{
			//printf("%d %d: %x %x %x %x\n", i, j, lambda_root, lambda[j], j * power_polynomial_table[i + 1][0], gf_multp(lambda[j], j * power_polynomial_table[i + 1][0]));
			lambda_root = gf_add(lambda_root, gf_multp(lambda[j], j * power_polynomial_table[i + 1][0]));
		}
		//printf("%d: %x %x\n", i, lambda_root, power_polynomial_table[i + 1][0]);
		if(0xFF == lambda_root)
		{
			err_location[GF_FIELD - 1 - i] = 1;/*inversion*/
		}
	}
	printf("Error Location:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", err_location[i]);
	}
	printf("\n");

	/*derivative of lambda*/
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2; i++)
	{
		if(0 != ((i + 1) % 2))
		{
			lambda_dev[i] = lambda[i + 1];
		}
	}
	printf("Lambda Derivative Polynomial:\n");
	for(i = 0; i < (CODEWORD_LEN - MESSAGE_LEN) / 2; i++)
	{
		printf("%x ", lambda_dev[i]);
	}
	printf("\n");

	/*magnitude of error*/
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		if(0 == err_location[i])
		{
			continue;
		}
		else
		{
			tmp = 0xFF;
			for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2 + (CODEWORD_LEN - MESSAGE_LEN - 1) / 2 + 1; j++)
			{
				//printf("%d %d: %x %x %x %x\n", i, j, tmp, omega[j], j * power_polynomial_table[i + 1][0], gf_multp(omega[j], j * power_polynomial_table[i + 1][0]));
				tmp = gf_add(tmp, gf_div(omega[j], j * power_polynomial_table[i + 1][0]));
			}
			tmp_sum = 0xFF;
			for(j = 0; j < (CODEWORD_LEN - MESSAGE_LEN) / 2; j++)
			{
				//printf("%d %d: %x %x %x %x\n", i, j, tmp_sum, lambda_dev[j], j * power_polynomial_table[i + 1][0], gf_multp(lambda_dev[j], j * power_polynomial_table[i + 1][0]));
				tmp_sum = gf_add(tmp_sum, gf_div(lambda_dev[j], j * power_polynomial_table[i + 1][0]));
			}
			//printf("%d: %x %x\n", i, tmp, tmp_sum);
			err_mag[i] = gf_div(tmp, tmp_sum);
		}
	}
	printf("Error Magnitude:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", err_mag[i]);
	}
	printf("\n");

	for(i = 0; i < CODEWORD_LEN; i++)
	{
		codeword[i] = gf_add(codeword[i], err_mag[i]);
	}
	printf("Codeword:\n");
	for(i = 0; i < CODEWORD_LEN; i++)
	{
		printf("%x ", codeword[i]);
	}
	printf("\n");

	return 0;
}
