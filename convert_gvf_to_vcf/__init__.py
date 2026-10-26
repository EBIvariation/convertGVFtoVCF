import os

env_memory_value = os.environ.get('MEMORY_DIAGNOSTIC')
if env_memory_value:
    MEMORY_DIAGNOSTIC = env_memory_value.strip().lower() in ('1', 'true', 'yes', 'on')
else:
    MEMORY_DIAGNOSTIC = False