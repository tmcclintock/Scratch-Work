section .data
	userMsg db 'Please enter a number: '
	lenUserMsg equ $-userMsg
	dispMsg db 'You have entered: '
	lenDispMsg equ $-dispMsg

section .bss
	num resb 5

section .text
	global _start

_start:
	mov eax,4		;sys_write
	mov ebx,1		;stdout
	mov ecx,userMsg
	mov edx,lenUserMsg
	int 0x80

	mov eax,3
	mov ebx,2
	mov ecx,num
	mov edx,5		;5 bytes (numeric, 1 for sign) of the num variable
	int 0x80
	
	mov eax,4		;sys_write
	mov ebx,1		;stdout
	mov ecx,dispMsg
	mov edx,lenDispMsg
	int 0x80

	mov eax,4
	mov ebx,1
	mov ecx,num
	mov edx,5
	int 0x80
	
	mov eax,1		;sys_exit
	int 0x80
