"""
Simple test script to demonstrate the comparison between PythonAstREPLTool and structured output approaches.
This is a minimal example that can be run quickly to test the basic functionality.
"""

import os
import pandas as pd
import scanpy as sc
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
    project_path = "test_project"
    project = MDVProject(project_path, delete_existing=True)
    project.add_datasource("test_data", df)
    
    return project

def test_pandas_agent(project: MDVProject, question: str):
    """Test the pandas agent approach"""
    print(f"\nTesting Pandas Agent: {question}")
    
    try:
        chat = ProjectChat(project)
        
        # Create chat request
        chat_request = ChatRequest(
            id="test_id",
            room="test_room", 
            conversation_id="test_conversation",
            message=question,
            handle_error=lambda error, **kwargs: print(f"Error: {error}")
        )
        
        # Run the query
        result = chat.ask_question(chat_request)
        
        print(f"Success: {not result.get('error', True)}")
        print(f"View name: {result.get('view_name', 'N/A')}")
        if result.get('error'):
            print(f"Error: {result.get('message', 'Unknown error')}")
        
        return result
        
    except Exception as e:
        print(f"Pandas Agent Error: {e}")
        return {"error": True, "message": str(e)}

def test_structured_approach(project: MDVProject, question: str):
    """Test the structured approach"""
    print(f"\nTesting Structured Approach: {question}")
    
    try:
        chat = StructuredProjectChat(project)
        
        # Create chat request
        chat_request = ChatRequest(
            id="test_id",
            room="test_room",
            conversation_id="test_conversation", 
            message=question,
            handle_error=lambda error, **kwargs: print(f"Error: {error}")
        )
        
        # Run the query
        result = chat.ask_question(chat_request)
        
        print(f"Success: {not result.get('error', True)}")
        print(f"View name: {result.get('view_name', 'N/A')}")
        if result.get('error'):
            print(f"Error: {result.get('message', 'Unknown error')}")
        
        # Print structured result if available
        if 'structured_result' in result:
            structured = result['structured_result']
            print(f"Data Analysis: {structured.get('data_analysis', {}).get('selected_columns', [])}")
            print(f"Chart Type: {structured.get('chart_generation', {}).get('chart_config', {}).get('type', 'N/A')}")
        
        return result
        
    except Exception as e:
        print(f"Structured Approach Error: {e}")
        return {"error": True, "message": str(e)}

def compare_approaches(project: MDVProject, questions: list):
    """Compare both approaches on a set of questions"""
    print("="*60)
    print("COMPARING PANDAS AGENT vs STRUCTURED APPROACH")
    print("="*60)
    
    results = []
    
    for i, question in enumerate(questions):
        print(f"\nQuestion {i+1}: {question}")
        print("-" * 40)
        
        # Test pandas agent
        pandas_result = test_pandas_agent(project, question)
        
        # Test structured approach
        structured_result = test_structured_approach(project, question)
        
        # Compare results
        pandas_success = not pandas_result.get('error', True)
        structured_success = not structured_result.get('error', True)
        
        print(f"\nComparison:")
        print(f"  Pandas Agent Success: {pandas_success}")
        print(f"  Structured Success: {structured_success}")
        print(f"  Success Match: {pandas_success == structured_success}")
        
        results.append({
            'question': question,
            'pandas_success': pandas_success,
            'structured_success': structured_success,
            'success_match': pandas_success == structured_success
        })
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    total_questions = len(results)
    pandas_successes = sum(1 for r in results if r['pandas_success'])
    structured_successes = sum(1 for r in results if r['structured_success'])
    success_matches = sum(1 for r in results if r['success_match'])
    
    print(f"Total Questions: {total_questions}")
    print(f"Pandas Agent Success Rate: {pandas_successes}/{total_questions} ({pandas_successes/total_questions:.1%})")
    print(f"Structured Success Rate: {structured_successes}/{total_questions} ({structured_successes/total_questions:.1%})")
    print(f"Success Match Rate: {success_matches}/{total_questions} ({success_matches/total_questions:.1%})")
    
    return results

def main():
    """Main function"""
    print("Simple LLM Approach Comparison Test")
    print("This test compares PythonAstREPLTool vs Structured Output approaches")
    
    # Check for OpenAI API key
    if not os.getenv("OPENAI_API_KEY"):
        print("ERROR: OPENAI_API_KEY not found in environment variables")
        print("Please set your OpenAI API key before running this test")
        return
    
    # Create test project
    project = create_test_project()
    
    # Test questions
    test_questions = [
        "Create a scatter plot of x vs y",
        "Show me a histogram of the value column",
        "Make a bar chart of categories"
    ]
    
    # Run comparison
    results = compare_approaches(project, test_questions)
    
    print("\nTest completed!")
    print("Check the test_project directory for generated visualizations")

if __name__ == "__main__":
    main() 