section	.text
   global _start     ;must be declared for linker (ld)
	
_start:	            ;tells linker entry point
   mov	edx,len     ;message length
   mov	ecx,msg     ;message to write
   mov	ebx,1       ;file descriptor (stdout)
   mov	eax,4       ;system call number (sys_write)
   int	0x80        ;call kernel

   mov  edx,9       ;message length
   mov  ecx,s2	    ;message to write
   mov  ebx,1	    ;file descriptor (stdout)
   mov  eax,4	    ;system call number (sys_write)
   int  0x80	    ;call kernel

   mov	eax,1       ;system call number (sys_exit)
   int	0x80        ;call kernel

section	.data
msg db 'Displaying 9 stars', 0xa  ;string to be printed
len equ $ - msg     ;length of the string
s2 times 9 db '*'
