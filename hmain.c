void __multchecks__(char* nbmults);


int main(int argc, char** argv)
{
	if(argc >= 2)
		__multchecks__(argv[1]);
	else
		__multchecks__('\0');
	return 0;
}
