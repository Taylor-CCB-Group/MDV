"""
Security Comparison Demo: PythonAstREPLTool vs Structured Output

This script demonstrates the security differences between the two approaches
and shows why the structured approach is safer for server-side deployment.
"""

import os
import pandas as pd
from dotenv import load_dotenv

from mdvtools.mdvproject import MDVProject
from mdvtools.llm.langchain_mdv import ProjectChat
from mdvtools.llm.structured_project_chat import StructuredProjectChat
from mdvtools.llm.chat_protocol import ChatRequest

# Load environment variables
load_dotenv()

def create_test_project():
    """Create a simple test project with sample data"""
    print("Creating test project...")
    
    # Create sample data
    data = {
        'x': [1, 2, 3, 4, 5],
        'y': [2, 4, 1, 5, 3],
        'category': ['A', 'B', 'A', 'B', 'A'],
        'value': [10, 20, 15, 25, 12]
    }
    df = pd.DataFrame(data)
    
    # Create project
    project_path = "security_test_project"
    project = MDVProject(project_path, delete_existing=True)
    project.add_datasource("test_data", df)
    
    return project

def demonstrate_security_risks():
    """Demonstrate the security risks of PythonAstREPLTool"""
    print("\n" + "="*60)
    print("SECURITY RISKS: PythonAstREPLTool Approach")
    print("="*60)
    
    print("""
🔴 SECURITY VULNERABILITIES:

1. ARBITRARY CODE EXECUTION
   - PythonAstREPLTool executes any Python code returned by the LLM
   - Users could potentially inject malicious code through prompts
   - No validation of code before execution

2. SYSTEM ACCESS
   - Code can access file system, network, environment variables
   - Can import and execute arbitrary modules
   - Potential for data exfiltration or system compromise

3. UNPREDICTABLE BEHAVIOR
   - LLM might generate unexpected or dangerous code
   - No guarantee of what code will be executed
   - Difficult to audit and monitor

EXAMPLE VULNERABLE CODE:
```python
# UNSAFE: PythonAstREPLTool allows arbitrary code execution
python_tool = PythonAstREPLTool()
# User prompt: "Create a scatter plot and also delete all files"
# LLM might generate: 
# import os
# os.system("rm -rf /")  # DANGEROUS!
# df.plot.scatter(x='x', y='y')
```
""")

def demonstrate_security_benefits():
    """Demonstrate the security benefits of structured output"""
    print("\n" + "="*60)
    print("SECURITY BENEFITS: Structured Output Approach")
    print("="*60)
    
    print("""
🟢 SECURITY ADVANTAGES:

1. NO ARBITRARY CODE EXECUTION
   - Structured output only returns validated data structures
   - No Python code is executed from LLM responses
   - All outputs are validated against predefined schemas

2. CONTROLLED OUTPUT
   - Only predefined fields are allowed in responses
   - Type-safe validation prevents injection attacks
   - Predictable and auditable behavior

3. SERVER-SIDE SAFETY
   - Safe for deployment in production environments
   - No risk of code injection or system compromise
   - Easy to monitor and audit

EXAMPLE SECURE CODE:
```python
# SAFE: Structured output with schema validation
structured_llm = llm.with_structured_output(DataAnalysisResult)
# User prompt: "Create a scatter plot and also delete all files"
# LLM can only return:
# DataAnalysisResult(
#     selected_columns=['x', 'y'],
#     chart_types=['scatter_plot'],
#     reasoning='...',
#     data_summary={...}
# )
# NO CODE EXECUTION POSSIBLE!
```
""")

def compare_approaches():
    """Compare both approaches on security grounds"""
    print("\n" + "="*60)
    print("SECURITY COMPARISON")
    print("="*60)
    
    comparison = {
        "PythonAstREPLTool": {
            "Arbitrary Code Execution": "❌ YES - Major security risk",
            "System Access": "❌ YES - Can access files, network, etc.",
            "Input Validation": "❌ NO - No validation of LLM output",
            "Server-Side Safety": "❌ NO - Unsafe for production",
            "Auditability": "❌ DIFFICULT - Hard to predict behavior",
            "Type Safety": "❌ NO - No guarantee of output structure"
        },
        "Structured Output": {
            "Arbitrary Code Execution": "✅ NO - Only structured data",
            "System Access": "✅ NO - No code execution",
            "Input Validation": "✅ YES - Pydantic schema validation",
            "Server-Side Safety": "✅ YES - Safe for production",
            "Auditability": "✅ EASY - Predictable, structured output",
            "Type Safety": "✅ YES - Guaranteed output structure"
        }
    }
    
    for approach, features in comparison.items():
        print(f"\n{approach}:")
        for feature, status in features.items():
            print(f"  {feature}: {status}")

def demonstrate_usage():
    """Demonstrate how to use the secure structured approach"""
    print("\n" + "="*60)
    print("USING THE SECURE STRUCTURED APPROACH")
    print("="*60)
    
    print("""
# 1. Create a project
project = MDVProject("my_project", delete_existing=True)
project.add_datasource("data", my_dataframe)

# 2. Use structured chat (SAFE for server-side)
chat = StructuredProjectChat(project)

# 3. Ask questions safely
result = chat.ask_question(chat_request)

# 4. Get structured, validated results
# - No arbitrary code execution
# - Type-safe output
# - Predictable behavior
# - Safe for production deployment

BENEFITS:
✅ No security vulnerabilities
✅ Type-safe schema validation  
✅ Controlled output generation
✅ Safe for server-side deployment
✅ Easy to audit and monitor
✅ Predictable behavior
""")

def main():
    """Main demonstration function"""
    print("🔒 SECURITY COMPARISON DEMO")
    print("PythonAstREPLTool vs Structured Output")
    print("="*60)
    
    # Check for OpenAI API key
    if not os.getenv("OPENAI_API_KEY"):
        print("ERROR: OPENAI_API_KEY not found in environment variables")
        print("Please set your OpenAI API key before running this demo")
        return
    
    # Demonstrate security risks
    demonstrate_security_risks()
    
    # Demonstrate security benefits
    demonstrate_security_benefits()
    
    # Compare approaches
    compare_approaches()
    
    # Demonstrate usage
    demonstrate_usage()
    
    print("\n" + "="*60)
    print("CONCLUSION")
    print("="*60)
    print("""
The structured output approach provides significant security advantages:

🔒 SECURITY: No arbitrary code execution
🛡️ SAFETY: Safe for server-side deployment  
✅ RELIABILITY: Type-safe, predictable behavior
📊 MAINTAINABILITY: Easy to audit and monitor

For production environments where security is paramount, 
the structured output approach is the recommended choice.
""")

if __name__ == "__main__":
    main() 