section .text
	global _start

_start:
	mov eax, 4		;sys_write
	mov ebx, 1		;stdout
	mov ecx, val
	mov edx, 1		;length
	int 0x80
	
	mov eax, 1
	int 0x80

section .data
val DD 'y'
