#include <stdio.h>
#include <gmp.h>

#define MODULUS_SIZE 2048
#define BLOCK_SIZE (MODULUS_SIZE/8)      
#define BUFFER_SIZE ((MODULUS_SIZE/8) / 2)
#define INPUT_FILE "input1.txt"

typedef struct {
	mpz_t n; /* Modulus */
	mpz_t e; /* Public Exponent */
} public_key;

typedef struct {
	mpz_t n; /* Modulus */
	mpz_t e; /* Public Exponent */
	mpz_t d; /* Private Exponent */
	mpz_t p; /* Starting prime p */
	mpz_t q; /* Starting prime q */
} private_key;

void print_hex(char* arr, int len)
{
	int i;
	for (i = 0; i < len; i++)
		printf("%02x", (unsigned char)arr[i]);
}

void generate_keys(private_key* ku, public_key* kp)
{
	char buf[BUFFER_SIZE];
	int i;
	mpz_t phi; mpz_init(phi);
	mpz_t tmp1; mpz_init(tmp1);
	mpz_t tmp2; mpz_init(tmp2);

	srand(time(NULL));

	/* Insetead of selecting e st. gcd(phi, e) = 1; 1 < e < phi, lets choose e
	* first then pick p,q st. gcd(e, p-1) = gcd(e, q-1) = 1 */
	// We'll set e globally.  I've seen suggestions to use primes like 3, 17 or 
	// 65537, as they make coming calculations faster.  Lets use 3.
	mpz_set_ui(ku->e, 65537);

	/* Select p and q */
	/* Start with p */
	// Set the bits of tmp randomly
	for (i = 0; i < BUFFER_SIZE; i++)
		buf[i] = rand() % 0xFF;
	// Set the top two bits to 1 to ensure int(tmp) is relatively large
	buf[0] |= 0xC0;
	// Set the bottom bit to 1 to ensure int(tmp) is odd (better for finding primes)
	buf[BUFFER_SIZE - 1] |= 0x01;
	// Interpret this char buffer as an int
	mpz_import(tmp1, BUFFER_SIZE, 1, sizeof(buf[0]), 0, 0, buf);
	// Pick the next prime starting from that random number
	mpz_nextprime(ku->p, tmp1);
	/* Make sure this is a good choice*/
	mpz_mod(tmp2, ku->p, ku->e);        /* If p mod e == 1, gcd(phi, e) != 1 */
	while (!mpz_cmp_ui(tmp2, 1))
	{
		mpz_nextprime(ku->p, ku->p);    /* so choose the next prime */
		mpz_mod(tmp2, ku->p, ku->e);
	}

	/* Now select q */
	do {
		for (i = 0; i < BUFFER_SIZE; i++)
			buf[i] = rand() % 0xFF;
		// Set the top two bits to 1 to ensure int(tmp) is relatively large
		buf[0] |= 0xC0;
		// Set the bottom bit to 1 to ensure int(tmp) is odd
		buf[BUFFER_SIZE - 1] |= 0x01;
		// Interpret this char buffer as an int
		mpz_import(tmp1, (BUFFER_SIZE), 1, sizeof(buf[0]), 0, 0, buf);
		// Pick the next prime starting from that random number
		mpz_nextprime(ku->q, tmp1);
		mpz_mod(tmp2, ku->q, ku->e);
		while (!mpz_cmp_ui(tmp2, 1))
		{
			mpz_nextprime(ku->q, ku->q);
			mpz_mod(tmp2, ku->q, ku->e);
		}
	} while (mpz_cmp(ku->p, ku->q) == 0); /* If we have identical primes (unlikely), try again */

										  /* Calculate n = p x q */
	mpz_mul(ku->n, ku->p, ku->q);

	/* Compute phi(n) = (p-1)(q-1) */
	mpz_sub_ui(tmp1, ku->p, 1);
	mpz_sub_ui(tmp2, ku->q, 1);
	mpz_mul(phi, tmp1, tmp2);

	/* Calculate d (multiplicative inverse of e mod phi) */
	if (mpz_invert(ku->d, ku->e, phi) == 0)
	{
		mpz_gcd(tmp1, ku->e, phi);
		printf("gcd(e, phi) = [%s]\n", mpz_get_str(NULL, 16, tmp1));
		printf("Invert failed\n");
	}

	/* Set public key */
	mpz_set(kp->e, ku->e);
	mpz_set(kp->n, ku->n);

	return;
}

void block_encrypt(mpz_t C, mpz_t M, public_key kp)
{
	/* C = M^e mod n */
	mpz_powm(C, M, kp.e, kp.n);
	return;
}

void block_decrypt(mpz_t M, mpz_t C, private_key ku)
{
	mpz_powm(M, C, ku.d, ku.n);
	return;
}

char *inputString(long long *numBlocks) {
	FILE *input;
	char c, *temp, *message;
	input = fopen(INPUT_FILE, "r");
	message = malloc(sizeof(char) * BLOCK_SIZE);
	if (!message) {
		printf("Error: Heap allocation failed.\n");
		return NULL;
	}
	int len = 0, blockidx = 0, numChars = 0;
	while (fscanf(input, "%c", &c) != EOF) {
		numChars++;
		*(message + blockidx * BLOCK_SIZE + len++) = c;
		if (len == BLOCK_SIZE) {
			temp = realloc(message, sizeof(char) * BLOCK_SIZE * (++blockidx + 1));
			if (temp) {
				message = temp;
				len = 0;
			}
			else {
				printf("Message reallocation failed!\n");
				return NULL;
			}
		}
	}
	fclose(input);
	*(message + blockidx * BLOCK_SIZE + len++) = '\0';
	temp = realloc(message, sizeof(char) * (BLOCK_SIZE * (blockidx + 1) + len));
	if (temp) {
		message = temp;
	}
	else {
		printf("Message reallocation failed!\n");
		return NULL;
	}
	*numBlocks = blockidx + 1;
	printf("Number of chars read: %d\n", numChars);
	return message;
}

int main() {
	int i;
	mpz_t M;  mpz_init(M);
	mpz_t C;  mpz_init(C);
	mpz_t DC;  mpz_init(DC);
	private_key ku;
	public_key kp;

	// Initialize public key
	mpz_init(kp.n);
	mpz_init(kp.e);
	// Initialize private key
	mpz_init(ku.n);
	mpz_init(ku.e);
	mpz_init(ku.d);
	mpz_init(ku.p);
	mpz_init(ku.q);

	generate_keys(&ku, &kp);
	printf("---------------Public Key-----------------\n");
	printf("kp.n is [%s]\n", mpz_get_str(NULL, 10, kp.n));
	printf("kp.e is [%s]\n", mpz_get_str(NULL, 10, kp.e));
	printf("---------------Private Key------------------\n");
	printf("ku.n is [%s]\n", mpz_get_str(NULL, 10, ku.n));
	printf("ku.e is [%s]\n", mpz_get_str(NULL, 10, ku.e));
	printf("ku.d is [%s]\n", mpz_get_str(NULL, 10, ku.d));
	printf("ku.p is [%s]\n", mpz_get_str(NULL, 10, ku.p));
	printf("ku.q is [%s]\n\n", mpz_get_str(NULL, 10, ku.q));

	char *message;
	int *numBlocks;
	numBlocks = malloc(sizeof(long long));
	message = inputString(numBlocks);
	if (message == NULL) {
		printf("Message read failed!\n");
		exit(1);
	}
	
	mpz_import(M, BLOCK_SIZE, 1, sizeof(char), 0, 0, msg);
	//mpz_import(M, (6 * BLOCK_SIZE), 1, sizeof(buf[0]), 0, 0, buf);
	printf("original is [%s]\n", mpz_get_str(NULL, 16, M));
	block_encrypt(C, M, kp);
	printf("encrypted is [%s]\n", mpz_get_str(NULL, 16, C));
	block_decrypt(DC, C, ku);
	printf("decrypted is [%s]\n", mpz_get_str(NULL, 16, DC));
	char result[BLOCK_SIZE+1];
	int length;
	mpz_export(result, &length, 1, sizeof(char), 0, 0, DC);
	result[BLOCK_SIZE] = '\0';
	printf("%s\nLength is %d\n", result, length);

	//getch();
	return 0;
}