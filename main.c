void __multchecks__(char* nbmults);
//void __sqandmultdemo(void);
//void __main__(void);

int main(int argc, char** argv)
{
	if(argc >= 2)
		__multchecks__(argv[1]);
	else
		__multchecks__('\0');
	return 0;
}
