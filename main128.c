void __multchecks__(char* nbmults);
void __sqandmultdemo(void);

int main(int argc, char** argv)
{
	if(argc >= 2)
		__multchecks__(argv[1]);
	else
		__multchecks__('\0');
	//__sqandmultdemo();
	return 0;
}
